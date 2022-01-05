/*
The compositions in the x() array use the following order and must be sent as mole fractions:
    0 - PLACEHOLDER

    1 - Methane
    2 - Nitrogen
    3 - Carbon dioxide
    4 - Ethane
    5 - Propane
    6 - Isobutane
    7 - n-Butane
    8 - Isopentane
    9 - n-Pentane
   10 - n-Hexane
   11 - n-Heptane
   12 - n-Octane
   13 - n-Nonane
   14 - n-Decane
   15 - Hydrogen
   16 - Oxygen
   17 - Carbon monoxide
   18 - Water
   19 - Hydrogen sulfide
   20 - Helium
   21 - Argon
*/
package main

import (
	"encoding/base64"
	"math"
	"strings"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/app"
	"fyne.io/fyne/v2/canvas"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/data/binding"
	"fyne.io/fyne/v2/layout"
	"fyne.io/fyne/v2/widget"
)

var RDetail float64

const (
	NcDetail = 21
	MaxFlds  = 21
	NTerms   = 58
)
const epsilon float64 = 1e-15

var (
	fn [NTerms + 1]int
	gn [NTerms + 1]int
	qn [NTerms + 1]int
)
var (
	an [NTerms + 1]float64
	un [NTerms + 1]float64
)
var (
	bn [NTerms + 1]int
	kn [NTerms + 1]int
)
var (
	Bsnij2 [MaxFlds + 1][MaxFlds + 1][18 + 1]float64
	Bs     [18 + 1]float64
	Csn    [NTerms + 1]float64
)
var (
	Fi [MaxFlds + 1]float64
	Gi [MaxFlds + 1]float64
	Qi [MaxFlds + 1]float64
)
var (
	Ki25 [MaxFlds + 1]float64
	Ei25 [MaxFlds + 1]float64
)
var (
	Kij5 [MaxFlds + 1][MaxFlds + 1]float64
	Uij5 [MaxFlds + 1][MaxFlds + 1]float64
	Gij5 [MaxFlds + 1][MaxFlds + 1]float64
)
var (
	Tun  [NTerms + 1]float64
	Told float64
)
var (
	n0i  [MaxFlds + 1][7 + 1]float64
	th0i [MaxFlds + 1][7 + 1]float64
)
var (
	MMiDetail [MaxFlds + 1]float64
	K3        float64
	xold      [MaxFlds + 1]float64
)
var dPdDsave float64

func sq(x float64) float64 { return x * x }

func MolarMassDetail(x []float64, Mm *float64) {
	// Calculate molar mass of the mixture with the compositions contained in the x() input array

	// Inputs:
	//    x() - Composition (mole fraction)
	//          Do not send mole percents or mass fractions in the x() array, otherwise the output will be incorrect.
	//          The sum of the compositions in the x() array must be equal to one.
	//          The order of the fluids in this array is given at the top of this code.

	// Outputs:
	//     Mm - Molar mass (g/mol)

	*Mm = 0
	for i := 1; i <= NcDetail; i++ {
		*Mm += x[i] * MMiDetail[i]
	}
}

func PressureDetail(T float64, D float64, x []float64, P *float64, Z *float64) {
	// Sub Pressure(T, D, x, P, Z)

	// Calculate pressure as a function of temperature and density.  The derivative d(P)/d(D) is also calculated
	// for use in the iterative DensityDetail subroutine (and is only returned as a common variable).

	// Inputs:
	//      T - Temperature (K)
	//      D - Density (mol/l)
	//    x() - Composition (mole fraction)
	//          Do not send mole percents or mass fractions in the x() array, otherwise the output will be incorrect.
	//          The sum of the compositions in the x() array must be equal to one.

	// Outputs:
	//      P - Pressure (kPa)
	//      Z - Compressibility factor
	//   dPdDsave - d(P)/d(D) [kPa/(mol/l)] (at constant temperature)
	//            - This variable is cached in the common variables for use in the iterative density solver, but not returned as an argument.

	var ar [3 + 1][3 + 1]float64
	xTermsDetail(x)
	AlpharDetail(0, 2, T, D, &ar)
	*Z = 1 + ar[0][1]/RDetail/T // ar(0,1) is the first derivative of alpha(r) with respect to density
	*P = D * RDetail * T * *Z
	dPdDsave = RDetail*T + 2*ar[0][1] + ar[0][2] // d(P)/d(D) for use in density iteration
}

func DensityDetail(T float64, P float64, x []float64, D *float64, ierr *int, herr *string) {
	// Sub DensityDetail(T, P, x, D, ierr, herr)

	// Calculate density as a function of temperature and pressure.  This is an iterative routine that calls PressureDetail
	// to find the correct state point.  Generally only 6 iterations at most are required.
	// If the iteration fails to converge, the ideal gas density and an error message are returned.
	// No checks are made to determine the phase boundary, which would have guaranteed that the output is in the gas phase.
	// It is up to the user to locate the phase boundary, and thus identify the phase of the T and P inputs.
	// If the state point is 2-phase, the output density will represent a metastable state.

	// Inputs:
	//      T - Temperature (K)
	//      P - Pressure (kPa)
	//    x() - Composition (mole fraction)

	// Outputs:
	//      D - Density (mol/l) (make D negative and send as an input to use as an initial guess)
	//   ierr - Error number (0 indicates no error)
	//   herr - Error message if ierr is not equal to zero

	var plog, vlog, P2, Z, dpdlv, vdiff, tolr float64

	*ierr = 0
	*herr = ""
	if math.Abs(P) < epsilon {
		*D = 0
		return
	}
	tolr = 0.0000001
	if *D > -epsilon {
		*D = P / RDetail / T // Ideal gas estimate
	} else {
		*D = math.Abs(*D) // If D<0, then use as initial estimate
	}
	plog = math.Log(P)
	vlog = -math.Log(*D)
	for it := 1; it <= 20; it++ {
		if vlog < -7 || vlog > 100 {
			*ierr = 1
			*herr = "Calculation failed to converge in DETAIL method, ideal gas density returned."
			*D = P / RDetail / T
			return
		}
		*D = math.Exp(-vlog)
		PressureDetail(T, *D, x, &P2, &Z)
		if dPdDsave < epsilon || P2 < epsilon {
			vlog += 0.1
		} else {
			// Find the next density with a first order Newton's type iterative scheme, with
			// log(P) as the known variable and log(v) as the unknown property.
			// See AGA 8 publication for further information.
			dpdlv = -*D * dPdDsave // d(p)/d[log(v)]
			vdiff = (math.Log(P2) - plog) * P2 / dpdlv
			vlog = vlog - vdiff
			if math.Abs(vdiff) < tolr {
				*D = math.Exp(-vlog)
				return // Iteration converged
			}
		}
	}
	*ierr = 1
	*herr = "Calculation failed to converge in DETAIL method, ideal gas density returned."
	*D = P / RDetail / T
}

func PropertiesDetail(T float64, D float64, x []float64, P *float64, Z *float64, dPdD *float64, d2PdD2 *float64, d2PdTD *float64, dPdT *float64, U *float64, H *float64, S *float64, Cv *float64, Cp *float64, W *float64, G *float64, JT *float64, Kappa *float64) {
	// Sub Properties(T, D, x, P, Z, dPdD, d2PdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa)

	// Calculate thermodynamic properties as a function of temperature and density.  Calls are made to the subroutines
	// Molarmass, Alpha0Detail, and AlpharDetail.  If the density is not known, call subroutine DensityDetail first
	// with the known values of pressure and temperature.

	// Inputs:
	//      T - Temperature (K)
	//      D - Density (mol/l)
	//    x() - Composition (mole fraction)

	// Outputs:
	//      P - Pressure (kPa)
	//      Z - Compressibility factor
	//   dPdD - First derivative of pressure with respect to density at constant temperature [kPa/(mol/l)]
	// d2PdD2 - Second derivative of pressure with respect to density at constant temperature [kPa/(mol/l)^2]
	// d2PdTD - Second derivative of pressure with respect to temperature and density [kPa/(mol/l)/K] (currently not calculated)
	//   dPdT - First derivative of pressure with respect to temperature at constant density (kPa/K)
	//      U - Internal energy (J/mol)
	//      H - Enthalpy (J/mol)
	//      S - Entropy [J/(mol-K)]
	//     Cv - Isochoric heat capacity [J/(mol-K)]
	//     Cp - Isobaric heat capacity [J/(mol-K)]
	//      W - Speed of sound (m/s)
	//      G - Gibbs energy (J/mol)
	//     JT - Joule-Thomson coefficient (K/kPa)
	//  Kappa - Isentropic Exponent

	var a0 [2 + 1]float64
	var ar [3 + 1][3 + 1]float64
	var Mm, A, R, RT float64

	MolarMassDetail(x, &Mm)
	xTermsDetail(x)

	// Calculate the ideal gas Helmholtz energy, and its first and second derivatives with respect to temperature.
	Alpha0Detail(T, D, x, &a0)

	// Calculate the real gas Helmholtz energy, and its derivatives with respect to temperature and/or density.
	AlpharDetail(2, 3, T, D, &ar)

	R = RDetail
	RT = R * T
	*Z = 1 + ar[0][1]/RT
	*P = D * RT * *Z
	*dPdD = RT + 2*ar[0][1] + ar[0][2]
	*dPdT = D*R + D*ar[1][1]
	A = a0[0] + ar[0][0]
	*S = -a0[1] - ar[1][0]
	*U = A + T**S
	*Cv = -(a0[2] + ar[2][0])
	if D > epsilon {
		*H = *U + *P/D
		*G = A + *P/D
		*Cp = *Cv + T*sq(*dPdT/D) / *dPdD
		*d2PdD2 = (2*ar[0][1] + 4*ar[0][2] + ar[0][3]) / D
		*JT = (T/D**dPdT / *dPdD - 1) / *Cp / D
	} else {
		*H = *U + RT
		*G = A + RT
		*Cp = *Cv + R
		*d2PdD2 = 0
		*JT = 1e+20 //=(dB/dT*T-B)/Cp for an ideal gas, but dB/dT is not calculated here
	}
	*W = 1000 * *Cp / *Cv * *dPdD / Mm
	if *W < 0 {
		*W = 0
	}
	*W = math.Sqrt(*W)
	*Kappa = *W * *W * Mm / (RT * 1000 * *Z)
	*d2PdTD = 0
}

// The following routines are low-level routines that should not be called outside of this code.
func xTermsDetail(x []float64) {
	// Calculate terms dependent only on composition
	//
	// Inputs:
	//    x() - Composition (mole fraction)

	var G, Q, F, U, Q2, xij, xi2 float64
	var icheck int

	// Check to see if a component fraction has changed.  If x is the same as the previous call, then exit.
	icheck = 0
	for i := 1; i <= NcDetail; i++ {
		if math.Abs(x[i]-xold[i]) > 0.0000001 {
			icheck = 1
		}
		xold[i] = x[i]
	}
	if icheck == 0 {
		return
	}

	K3 = 0
	U = 0
	G = 0
	Q = 0
	F = 0
	for n := 1; n <= 18; n++ {
		Bs[n] = 0
	}

	// Calculate pure fluid contributions
	for i := 1; i <= NcDetail; i++ {
		if x[i] > 0 {
			xi2 = sq(x[i])
			K3 += x[i] * Ki25[i] // K, U, and G are the sums of a pure fluid contribution and a
			U += x[i] * Ei25[i]  // binary pair contribution
			G += x[i] * Gi[i]
			Q += x[i] * Qi[i] // Q and F depend only on the pure fluid parts
			F += xi2 * Fi[i]
			for n := 1; n <= 18; n++ {
				Bs[n] = Bs[n] + xi2*Bsnij2[i][i][n] // Pure fluid contributions to second virial coefficient
			}
		}
	}
	K3 = sq(K3)
	U = sq(U)

	// Binary pair contributions
	for i := 1; i <= NcDetail-1; i++ {
		if x[i] > 0 {
			for j := i + 1; j <= NcDetail; j++ {
				if x[j] > 0 {
					xij = 2 * x[i] * x[j]
					K3 = K3 + xij*Kij5[i][j]
					U = U + xij*Uij5[i][j]
					G = G + xij*Gij5[i][j]
					for n := 1; n <= 18; n++ {
						Bs[n] = Bs[n] + xij*Bsnij2[i][j][n] // Second virial coefficients of mixture
					}
				}
			}
		}
	}
	K3 = math.Pow(K3, 0.6)
	U = math.Pow(U, 0.2)

	// Third virial and higher coefficients
	Q2 = sq(Q)
	for n := 13; n <= 58; n++ {
		Csn[n] = an[n] * math.Pow(U, un[n])
		if gn[n] == 1 {
			Csn[n] = Csn[n] * G
		}
		if qn[n] == 1 {
			Csn[n] = Csn[n] * Q2
		}
		if fn[n] == 1 {
			Csn[n] = Csn[n] * F
		}
	}
}

func Alpha0Detail(T float64, D float64, x []float64, a0 *[3]float64) {
	// Private Sub Alpha0Detail(T, D, x, a0)

	// Calculate the ideal gas Helmholtz energy and its derivatives with respect to T and D.
	// This routine is not needed when only P (or Z) is calculated.

	// Inputs:
	//      T - Temperature (K)
	//      D - Density (mol/l)
	//    x() - Composition (mole fraction)

	// Outputs:
	// a0(0) - Ideal gas Helmholtz energy (J/mol)
	// a0(1) -   partial  (a0)/partial(T) [J/(mol-K)]
	// a0(2) - T*partial^2(a0)/partial(T)^2 [J/(mol-K)]

	var LogT, LogD, LogHyp, th0T, LogxD float64
	var SumHyp0, SumHyp1, SumHyp2 float64
	var em, ep, hcn, hsn float64

	a0[0] = 0
	a0[1] = 0
	a0[2] = 0
	if D > epsilon {
		LogD = math.Log(D)
	} else {
		LogD = math.Log(epsilon)
	}
	LogT = math.Log(T)
	for i := 1; i <= NcDetail; i++ {
		if x[i] > 0 {
			LogxD = LogD + math.Log(x[i])
			SumHyp0 = 0
			SumHyp1 = 0
			SumHyp2 = 0
			for j := 4; j <= 7; j++ {
				if th0i[i][j] > 0 {
					th0T = th0i[i][j] / T
					ep = math.Exp(th0T)
					em = 1 / ep
					hsn = (ep - em) / 2
					hcn = (ep + em) / 2
					if j == 4 || j == 6 {
						LogHyp = math.Log(math.Abs(hsn))
						SumHyp0 += n0i[i][j] * LogHyp
						SumHyp1 += n0i[i][j] * (LogHyp - th0T*hcn/hsn)
						SumHyp2 += n0i[i][j] * sq(th0T/hsn)
					} else {
						LogHyp = math.Log(math.Abs(hcn))
						SumHyp0 += -n0i[i][j] * LogHyp
						SumHyp1 += -n0i[i][j] * (LogHyp - th0T*hsn/hcn)
						SumHyp2 += +n0i[i][j] * sq(th0T/hcn)
					}
				}
			}
			a0[0] += x[i] * (LogxD + n0i[i][1] + n0i[i][2]/T - n0i[i][3]*LogT + SumHyp0)
			a0[1] += x[i] * (LogxD + n0i[i][1] - n0i[i][3]*(1+LogT) + SumHyp1)
			a0[2] += -x[i] * (n0i[i][3] + SumHyp2)
		}
	}
	a0[0] = a0[0] * RDetail * T
	a0[1] = a0[1] * RDetail
	a0[2] = a0[2] * RDetail
}

func AlpharDetail(itau int, idel int, T float64, D float64, ar *[4][4]float64) {
	// Private Sub AlpharDetail(itau, idel, T, D, ar)

	// Calculate the derivatives of the residual Helmholtz energy (ar) with respect to T and D.
	// itau and idel are inputs that contain the highest derivatives needed.
	// Outputs are returned in the array ar.
	// Subroutine xTerms must be called before this routine if x has changed

	// Inputs:
	//  itau - Set this to 1 to calculate "ar" derivatives with respect to T [i.e., ar(1,0), ar(1,1), and ar(2,0)], otherwise set it to 0.
	//  idel - Currently not used, but kept as an input for future use in specifing the highest density derivative needed.
	//     T - Temperature (K)
	//     D - Density (mol/l)

	// Outputs:
	// ar(0,0) - Residual Helmholtz energy (J/mol)
	// ar(0,1) -   D*partial  (ar)/partial(D) (J/mol)
	// ar(0,2) - D^2*partial^2(ar)/partial(D)^2 (J/mol)
	// ar(0,3) - D^3*partial^3(ar)/partial(D)^3 (J/mol)
	// ar(1,0) -     partial  (ar)/partial(T) [J/(mol-K)]
	// ar(1,1) -   D*partial^2(ar)/partial(D)/partial(T) [J/(mol-K)]
	// ar(2,0) -   T*partial^2(ar)/partial(T)^2 [J/(mol-K)]

	var ckd, bkd, Dred float64
	var Sum, s0, s1, s2, s3, RT float64
	var (
		Sum0 [NTerms + 1]float64
		SumB [NTerms + 1]float64
		Dknn [9 + 1]float64
		Expn [4 + 1]float64
	)
	var (
		CoefD1 [NTerms + 1]float64
		CoefD2 [NTerms + 1]float64
		CoefD3 [NTerms + 1]float64
	)
	var (
		CoefT1 [NTerms + 1]float64
		CoefT2 [NTerms + 1]float64
	)

	for i := 0; i <= 3; i++ {
		for j := 0; j <= 3; j++ {
			ar[i][j] = 0
		}
	}
	if math.Abs(T-Told) > 0.0000001 {
		for n := 1; n <= 58; n++ {
			Tun[n] = math.Pow(T, -un[n])
		}
	}
	Told = T

	// Precalculation of common powers and exponents of density
	Dred = K3 * D
	Dknn[0] = 1
	for n := 1; n <= 9; n++ {
		Dknn[n] = Dred * Dknn[n-1]
	}
	Expn[0] = 1
	for n := 1; n <= 4; n++ {
		Expn[n] = math.Exp(-Dknn[n])
	}
	RT = RDetail * T

	for n := 1; n <= 58; n++ {
		// Contributions to the Helmholtz energy and its derivatives with respect to temperature
		CoefT1[n] = RDetail * (un[n] - 1)
		CoefT2[n] = CoefT1[n] * un[n]
		// Contributions to the virial coefficients
		SumB[n] = 0
		Sum0[n] = 0
		if n <= 18 {
			Sum = Bs[n] * D
			if n >= 13 {
				Sum += -Csn[n] * Dred
			}
			SumB[n] = Sum * Tun[n]
		}
		if n >= 13 {
			// Contributions to the residual part of the Helmholtz energy
			Sum0[n] = Csn[n] * Dknn[bn[n]] * Tun[n] * Expn[kn[n]]
			// Contributions to the derivatives of the Helmholtz energy with respect to density
			bkd = float64(bn[n]) - float64(kn[n])*Dknn[kn[n]]
			ckd = float64(kn[n]) * float64(kn[n]) * Dknn[kn[n]]
			CoefD1[n] = bkd
			CoefD2[n] = bkd*(bkd-1) - ckd
			CoefD3[n] = (bkd-2)*CoefD2[n] + ckd*(1-float64(kn[n])-2*bkd)
		} else {
			CoefD1[n] = 0
			CoefD2[n] = 0
			CoefD3[n] = 0
		}
	}

	for n := 1; n <= 58; n++ {
		// Density derivatives
		s0 = Sum0[n] + SumB[n]
		s1 = Sum0[n]*CoefD1[n] + SumB[n]
		s2 = Sum0[n] * CoefD2[n]
		s3 = Sum0[n] * CoefD3[n]
		ar[0][0] = ar[0][0] + RT*s0
		ar[0][1] = ar[0][1] + RT*s1
		ar[0][2] = ar[0][2] + RT*s2
		ar[0][3] = ar[0][3] + RT*s3
		// Temperature derivatives
		if itau > 0 {
			ar[1][0] = ar[1][0] - CoefT1[n]*s0
			ar[1][1] = ar[1][1] - CoefT1[n]*s1
			ar[2][0] = ar[2][0] + CoefT2[n]*s0
		}
	}
}

/// The following routine must be called once before any other routine.
func SetupDetail() {
	// Initialize all the constants and parameters in the DETAIL model.
	// Some values are modified for calculations that do not depend on T, D, and x in order to speed up the program.

	var sn [NTerms + 1]int
	var wn [NTerms + 1]int
	var (
		Ei    [MaxFlds + 1]float64
		Ki    [MaxFlds + 1]float64
		Si    [MaxFlds + 1]float64
		Wi    [MaxFlds + 1]float64
		Bsnij float64
	)
	var (
		Kij [MaxFlds + 1][MaxFlds + 1]float64
		Gij [MaxFlds + 1][MaxFlds + 1]float64
		Eij [MaxFlds + 1][MaxFlds + 1]float64
		Uij [MaxFlds + 1][MaxFlds + 1]float64
	)
	var d0 float64

	RDetail = 8.31451

	// Molar masses (g/mol)
	MMiDetail[1] = 16.043   // Methane
	MMiDetail[2] = 28.0135  // Nitrogen
	MMiDetail[3] = 44.01    // Carbon dioxide
	MMiDetail[4] = 30.07    // Ethane
	MMiDetail[5] = 44.097   // Propane
	MMiDetail[6] = 58.123   // Isobutane
	MMiDetail[7] = 58.123   // n-Butane
	MMiDetail[8] = 72.15    // Isopentane
	MMiDetail[9] = 72.15    // n-Pentane
	MMiDetail[10] = 86.177  // Hexane
	MMiDetail[11] = 100.204 // Heptane
	MMiDetail[12] = 114.231 // Octane
	MMiDetail[13] = 128.258 // Nonane
	MMiDetail[14] = 142.285 // Decane
	MMiDetail[15] = 2.0159  // Hydrogen
	MMiDetail[16] = 31.9988 // Oxygen
	MMiDetail[17] = 28.01   // Carbon monoxide
	MMiDetail[18] = 18.0153 // Water
	MMiDetail[19] = 34.082  // Hydrogen sulfide
	MMiDetail[20] = 4.0026  // Helium
	MMiDetail[21] = 39.948  // Argon

	// Initialize constants
	Told = 0
	for i := 1; i <= NTerms; i++ {
		an[i] = 0
		bn[i] = 0
		gn[i] = 0
		fn[i] = 0
		kn[i] = 0
		qn[i] = 0
		sn[i] = 0
		un[i] = 0
		wn[i] = 0
	}

	for i := 1; i <= MaxFlds; i++ {
		Ei[i] = 0
		Fi[i] = 0
		Gi[i] = 0
		Ki[i] = 0
		Qi[i] = 0
		Si[i] = 0
		Wi[i] = 0
		xold[i] = 0
		for j := 1; j <= MaxFlds; j++ {
			Eij[i][j] = 1
			Gij[i][j] = 1
			Kij[i][j] = 1
			Uij[i][j] = 1
		}
	}

	// Coefficients of the equation of state
	an[1] = 0.1538326
	an[2] = 1.341953
	an[3] = -2.998583
	an[4] = -0.04831228
	an[5] = 0.3757965
	an[6] = -1.589575
	an[7] = -0.05358847
	an[8] = 0.88659463
	an[9] = -0.71023704
	an[10] = -1.471722
	an[11] = 1.32185035
	an[12] = -0.78665925
	an[13] = 0.00000000229129
	an[14] = 0.1576724
	an[15] = -0.4363864
	an[16] = -0.04408159
	an[17] = -0.003433888
	an[18] = 0.03205905
	an[19] = 0.02487355
	an[20] = 0.07332279
	an[21] = -0.001600573
	an[22] = 0.6424706
	an[23] = -0.4162601
	an[24] = -0.06689957
	an[25] = 0.2791795
	an[26] = -0.6966051
	an[27] = -0.002860589
	an[28] = -0.008098836
	an[29] = 3.150547
	an[30] = 0.007224479
	an[31] = -0.7057529
	an[32] = 0.5349792
	an[33] = -0.07931491
	an[34] = -1.418465
	an[35] = -5.99905e-17
	an[36] = 0.1058402
	an[37] = 0.03431729
	an[38] = -0.007022847
	an[39] = 0.02495587
	an[40] = 0.04296818
	an[41] = 0.7465453
	an[42] = -0.2919613
	an[43] = 7.294616
	an[44] = -9.936757
	an[45] = -0.005399808
	an[46] = -0.2432567
	an[47] = 0.04987016
	an[48] = 0.003733797
	an[49] = 1.874951
	an[50] = 0.002168144
	an[51] = -0.6587164
	an[52] = 0.000205518
	an[53] = 0.009776195
	an[54] = -0.02048708
	an[55] = 0.01557322
	an[56] = 0.006862415
	an[57] = -0.001226752
	an[58] = 0.002850908

	// Density exponents
	bn[1] = 1
	bn[2] = 1
	bn[3] = 1
	bn[4] = 1
	bn[5] = 1
	bn[6] = 1
	bn[7] = 1
	bn[8] = 1
	bn[9] = 1
	bn[10] = 1
	bn[11] = 1
	bn[12] = 1
	bn[13] = 1
	bn[14] = 1
	bn[15] = 1
	bn[16] = 1
	bn[17] = 1
	bn[18] = 1
	bn[19] = 2
	bn[20] = 2
	bn[21] = 2
	bn[22] = 2
	bn[23] = 2
	bn[24] = 2
	bn[25] = 2
	bn[26] = 2
	bn[27] = 2
	bn[28] = 3
	bn[29] = 3
	bn[30] = 3
	bn[31] = 3
	bn[32] = 3
	bn[33] = 3
	bn[34] = 3
	bn[35] = 3
	bn[36] = 3
	bn[37] = 3
	bn[38] = 4
	bn[39] = 4
	bn[40] = 4
	bn[41] = 4
	bn[42] = 4
	bn[43] = 4
	bn[44] = 4
	bn[45] = 5
	bn[46] = 5
	bn[47] = 5
	bn[48] = 5
	bn[49] = 5
	bn[50] = 6
	bn[51] = 6
	bn[52] = 7
	bn[53] = 7
	bn[54] = 8
	bn[55] = 8
	bn[56] = 8
	bn[57] = 9
	bn[58] = 9

	// Exponents on density in EXP[-cn*D^kn] part
	// The cn part in this term is not included in this program since it is 1 when kn<>0][and 0 otherwise
	kn[13] = 3
	kn[14] = 2
	kn[15] = 2
	kn[16] = 2
	kn[17] = 4
	kn[18] = 4
	kn[21] = 2
	kn[22] = 2
	kn[23] = 2
	kn[24] = 4
	kn[25] = 4
	kn[26] = 4
	kn[27] = 4
	kn[29] = 1
	kn[30] = 1
	kn[31] = 2
	kn[32] = 2
	kn[33] = 3
	kn[34] = 3
	kn[35] = 4
	kn[36] = 4
	kn[37] = 4
	kn[40] = 2
	kn[41] = 2
	kn[42] = 2
	kn[43] = 4
	kn[44] = 4
	kn[46] = 2
	kn[47] = 2
	kn[48] = 4
	kn[49] = 4
	kn[51] = 2
	kn[53] = 2
	kn[54] = 1
	kn[55] = 2
	kn[56] = 2
	kn[57] = 2
	kn[58] = 2

	// Temperature exponents
	un[1] = 0
	un[2] = 0.5
	un[3] = 1
	un[4] = 3.5
	un[5] = -0.5
	un[6] = 4.5
	un[7] = 0.5
	un[8] = 7.5
	un[9] = 9.5
	un[10] = 6
	un[11] = 12
	un[12] = 12.5
	un[13] = -6
	un[14] = 2
	un[15] = 3
	un[16] = 2
	un[17] = 2
	un[18] = 11
	un[19] = -0.5
	un[20] = 0.5
	un[21] = 0
	un[22] = 4
	un[23] = 6
	un[24] = 21
	un[25] = 23
	un[26] = 22
	un[27] = -1
	un[28] = -0.5
	un[29] = 7
	un[30] = -1
	un[31] = 6
	un[32] = 4
	un[33] = 1
	un[34] = 9
	un[35] = -13
	un[36] = 21
	un[37] = 8
	un[38] = -0.5
	un[39] = 0
	un[40] = 2
	un[41] = 7
	un[42] = 9
	un[43] = 22
	un[44] = 23
	un[45] = 1
	un[46] = 9
	un[47] = 3
	un[48] = 8
	un[49] = 23
	un[50] = 1.5
	un[51] = 5
	un[52] = -0.5
	un[53] = 4
	un[54] = 7
	un[55] = 3
	un[56] = 0
	un[57] = 1
	un[58] = 0

	// Flags
	fn[13] = 1
	fn[27] = 1
	fn[30] = 1
	fn[35] = 1
	gn[5] = 1
	gn[6] = 1
	gn[25] = 1
	gn[29] = 1
	gn[32] = 1
	gn[33] = 1
	gn[34] = 1
	gn[51] = 1
	gn[54] = 1
	gn[56] = 1
	qn[7] = 1
	qn[16] = 1
	qn[26] = 1
	qn[28] = 1
	qn[37] = 1
	qn[42] = 1
	qn[47] = 1
	qn[49] = 1
	qn[52] = 1
	qn[58] = 1
	sn[8] = 1
	sn[9] = 1
	wn[10] = 1
	wn[11] = 1
	wn[12] = 1

	// Energy parameters
	Ei[1] = 151.3183
	Ei[2] = 99.73778
	Ei[3] = 241.9606
	Ei[4] = 244.1667
	Ei[5] = 298.1183
	Ei[6] = 324.0689
	Ei[7] = 337.6389
	Ei[8] = 365.5999
	Ei[9] = 370.6823
	Ei[10] = 402.636293
	Ei[11] = 427.72263
	Ei[12] = 450.325022
	Ei[13] = 470.840891
	Ei[14] = 489.558373
	Ei[15] = 26.95794
	Ei[16] = 122.7667
	Ei[17] = 105.5348
	Ei[18] = 514.0156
	Ei[19] = 296.355
	Ei[20] = 2.610111
	Ei[21] = 119.6299

	// Size parameters
	Ki[1] = 0.4619255
	Ki[2] = 0.4479153
	Ki[3] = 0.4557489
	Ki[4] = 0.5279209
	Ki[5] = 0.583749
	Ki[6] = 0.6406937
	Ki[7] = 0.6341423
	Ki[8] = 0.6738577
	Ki[9] = 0.6798307
	Ki[10] = 0.7175118
	Ki[11] = 0.7525189
	Ki[12] = 0.784955
	Ki[13] = 0.8152731
	Ki[14] = 0.8437826
	Ki[15] = 0.3514916
	Ki[16] = 0.4186954
	Ki[17] = 0.4533894
	Ki[18] = 0.3825868
	Ki[19] = 0.4618263
	Ki[20] = 0.3589888
	Ki[21] = 0.4216551

	// Orientation parameters
	Gi[2] = 0.027815
	Gi[3] = 0.189065
	Gi[4] = 0.0793
	Gi[5] = 0.141239
	Gi[6] = 0.256692
	Gi[7] = 0.281835
	Gi[8] = 0.332267
	Gi[9] = 0.366911
	Gi[10] = 0.289731
	Gi[11] = 0.337542
	Gi[12] = 0.383381
	Gi[13] = 0.427354
	Gi[14] = 0.469659
	Gi[15] = 0.034369
	Gi[16] = 0.021
	Gi[17] = 0.038953
	Gi[18] = 0.3325
	Gi[19] = 0.0885

	// Quadrupole parameters
	Qi[3] = 0.69
	Qi[18] = 1.06775
	Qi[19] = 0.633276
	Fi[15] = 1      // High temperature parameter
	Si[18] = 1.5822 // Dipole parameter
	Si[19] = 0.39   // Dipole parameter
	Wi[18] = 1      // Association parameter

	// Energy parameters
	Eij[1][2] = 0.97164
	Eij[1][3] = 0.960644
	Eij[1][5] = 0.994635
	Eij[1][6] = 1.01953
	Eij[1][7] = 0.989844
	Eij[1][8] = 1.00235
	Eij[1][9] = 0.999268
	Eij[1][10] = 1.107274
	Eij[1][11] = 0.88088
	Eij[1][12] = 0.880973
	Eij[1][13] = 0.881067
	Eij[1][14] = 0.881161
	Eij[1][15] = 1.17052
	Eij[1][17] = 0.990126
	Eij[1][18] = 0.708218
	Eij[1][19] = 0.931484
	Eij[2][3] = 1.02274
	Eij[2][4] = 0.97012
	Eij[2][5] = 0.945939
	Eij[2][6] = 0.946914
	Eij[2][7] = 0.973384
	Eij[2][8] = 0.95934
	Eij[2][9] = 0.94552
	Eij[2][15] = 1.08632
	Eij[2][16] = 1.021
	Eij[2][17] = 1.00571
	Eij[2][18] = 0.746954
	Eij[2][19] = 0.902271
	Eij[3][4] = 0.925053
	Eij[3][5] = 0.960237
	Eij[3][6] = 0.906849
	Eij[3][7] = 0.897362
	Eij[3][8] = 0.726255
	Eij[3][9] = 0.859764
	Eij[3][10] = 0.855134
	Eij[3][11] = 0.831229
	Eij[3][12] = 0.80831
	Eij[3][13] = 0.786323
	Eij[3][14] = 0.765171
	Eij[3][15] = 1.28179
	Eij[3][17] = 1.5
	Eij[3][18] = 0.849408
	Eij[3][19] = 0.955052
	Eij[4][5] = 1.02256
	Eij[4][7] = 1.01306
	Eij[4][9] = 1.00532
	Eij[4][15] = 1.16446
	Eij[4][18] = 0.693168
	Eij[4][19] = 0.946871
	Eij[5][7] = 1.0049
	Eij[5][15] = 1.034787
	Eij[6][15] = 1.3
	Eij[7][15] = 1.3
	Eij[10][19] = 1.008692
	Eij[11][19] = 1.010126
	Eij[12][19] = 1.011501
	Eij[13][19] = 1.012821
	Eij[14][19] = 1.014089
	Eij[15][17] = 1.1

	// Conformal energy parameters
	Uij[1][2] = 0.886106
	Uij[1][3] = 0.963827
	Uij[1][5] = 0.990877
	Uij[1][7] = 0.992291
	Uij[1][9] = 1.00367
	Uij[1][10] = 1.302576
	Uij[1][11] = 1.191904
	Uij[1][12] = 1.205769
	Uij[1][13] = 1.219634
	Uij[1][14] = 1.233498
	Uij[1][15] = 1.15639
	Uij[1][19] = 0.736833
	Uij[2][3] = 0.835058
	Uij[2][4] = 0.816431
	Uij[2][5] = 0.915502
	Uij[2][7] = 0.993556
	Uij[2][15] = 0.408838
	Uij[2][19] = 0.993476
	Uij[3][4] = 0.96987
	Uij[3][10] = 1.066638
	Uij[3][11] = 1.077634
	Uij[3][12] = 1.088178
	Uij[3][13] = 1.098291
	Uij[3][14] = 1.108021
	Uij[3][17] = 0.9
	Uij[3][19] = 1.04529
	Uij[4][5] = 1.065173
	Uij[4][6] = 1.25
	Uij[4][7] = 1.25
	Uij[4][8] = 1.25
	Uij[4][9] = 1.25
	Uij[4][15] = 1.61666
	Uij[4][19] = 0.971926
	Uij[10][19] = 1.028973
	Uij[11][19] = 1.033754
	Uij[12][19] = 1.038338
	Uij[13][19] = 1.042735
	Uij[14][19] = 1.046966

	// Size parameters
	Kij[1][2] = 1.00363
	Kij[1][3] = 0.995933
	Kij[1][5] = 1.007619
	Kij[1][7] = 0.997596
	Kij[1][9] = 1.002529
	Kij[1][10] = 0.982962
	Kij[1][11] = 0.983565
	Kij[1][12] = 0.982707
	Kij[1][13] = 0.981849
	Kij[1][14] = 0.980991
	Kij[1][15] = 1.02326
	Kij[1][19] = 1.00008
	Kij[2][3] = 0.982361
	Kij[2][4] = 1.00796
	Kij[2][15] = 1.03227
	Kij[2][19] = 0.942596
	Kij[3][4] = 1.00851
	Kij[3][10] = 0.910183
	Kij[3][11] = 0.895362
	Kij[3][12] = 0.881152
	Kij[3][13] = 0.86752
	Kij[3][14] = 0.854406
	Kij[3][19] = 1.00779
	Kij[4][5] = 0.986893
	Kij[4][15] = 1.02034
	Kij[4][19] = 0.999969
	Kij[10][19] = 0.96813
	Kij[11][19] = 0.96287
	Kij[12][19] = 0.957828
	Kij[13][19] = 0.952441
	Kij[14][19] = 0.948338

	// Orientation parameters
	Gij[1][3] = 0.807653
	Gij[1][15] = 1.95731
	Gij[2][3] = 0.982746
	Gij[3][4] = 0.370296
	Gij[3][18] = 1.67309

	// Ideal gas parameters
	n0i[1][3] = 4.00088
	n0i[1][4] = 0.76315
	n0i[1][5] = 0.0046
	n0i[1][6] = 8.74432
	n0i[1][7] = -4.46921
	n0i[1][1] = 29.83843397
	n0i[1][2] = -15999.69151
	n0i[2][3] = 3.50031
	n0i[2][4] = 0.13732
	n0i[2][5] = -0.1466
	n0i[2][6] = 0.90066
	n0i[2][7] = 0
	n0i[2][1] = 17.56770785
	n0i[2][2] = -2801.729072
	n0i[3][3] = 3.50002
	n0i[3][4] = 2.04452
	n0i[3][5] = -1.06044
	n0i[3][6] = 2.03366
	n0i[3][7] = 0.01393
	n0i[3][1] = 20.65844696
	n0i[3][2] = -4902.171516
	n0i[4][3] = 4.00263
	n0i[4][4] = 4.33939
	n0i[4][5] = 1.23722
	n0i[4][6] = 13.1974
	n0i[4][7] = -6.01989
	n0i[4][1] = 36.73005938
	n0i[4][2] = -23639.65301
	n0i[5][3] = 4.02939
	n0i[5][4] = 6.60569
	n0i[5][5] = 3.197
	n0i[5][6] = 19.1921
	n0i[5][7] = -8.37267
	n0i[5][1] = 44.70909619
	n0i[5][2] = -31236.63551
	n0i[6][3] = 4.06714
	n0i[6][4] = 8.97575
	n0i[6][5] = 5.25156
	n0i[6][6] = 25.1423
	n0i[6][7] = 16.1388
	n0i[6][1] = 34.30180349
	n0i[6][2] = -38525.50276
	n0i[7][3] = 4.33944
	n0i[7][4] = 9.44893
	n0i[7][5] = 6.89406
	n0i[7][6] = 24.4618
	n0i[7][7] = 14.7824
	n0i[7][1] = 36.53237783
	n0i[7][2] = -38957.80933
	n0i[8][3] = 4
	n0i[8][4] = 11.7618
	n0i[8][5] = 20.1101
	n0i[8][6] = 33.1688
	n0i[8][7] = 0
	n0i[8][1] = 43.17218626
	n0i[8][2] = -51198.30946
	n0i[9][3] = 4
	n0i[9][4] = 8.95043
	n0i[9][5] = 21.836
	n0i[9][6] = 33.4032
	n0i[9][7] = 0
	n0i[9][1] = 42.67837089
	n0i[9][2] = -45215.83
	n0i[10][3] = 4
	n0i[10][4] = 11.6977
	n0i[10][5] = 26.8142
	n0i[10][6] = 38.6164
	n0i[10][7] = 0
	n0i[10][1] = 46.99717188
	n0i[10][2] = -52746.83318
	n0i[11][3] = 4
	n0i[11][4] = 13.7266
	n0i[11][5] = 30.4707
	n0i[11][6] = 43.5561
	n0i[11][7] = 0
	n0i[11][1] = 52.07631631
	n0i[11][2] = -57104.81056
	n0i[12][3] = 4
	n0i[12][4] = 15.6865
	n0i[12][5] = 33.8029
	n0i[12][6] = 48.1731
	n0i[12][7] = 0
	n0i[12][1] = 57.25830934
	n0i[12][2] = -60546.76385
	n0i[13][3] = 4
	n0i[13][4] = 18.0241
	n0i[13][5] = 38.1235
	n0i[13][6] = 53.3415
	n0i[13][7] = 0
	n0i[13][1] = 62.09646901
	n0i[13][2] = -66600.12837
	n0i[14][3] = 4
	n0i[14][4] = 21.0069
	n0i[14][5] = 43.4931
	n0i[14][6] = 58.3657
	n0i[14][7] = 0
	n0i[14][1] = 65.93909154
	n0i[14][2] = -74131.45483
	n0i[15][3] = 2.47906
	n0i[15][4] = 0.95806
	n0i[15][5] = 0.45444
	n0i[15][6] = 1.56039
	n0i[15][7] = -1.3756
	n0i[15][1] = 13.07520288
	n0i[15][2] = -5836.943696
	n0i[16][3] = 3.50146
	n0i[16][4] = 1.07558
	n0i[16][5] = 1.01334
	n0i[16][6] = 0
	n0i[16][7] = 0
	n0i[16][1] = 16.8017173
	n0i[16][2] = -2318.32269
	n0i[17][3] = 3.50055
	n0i[17][4] = 1.02865
	n0i[17][5] = 0.00493
	n0i[17][6] = 0
	n0i[17][7] = 0
	n0i[17][1] = 17.45786899
	n0i[17][2] = -2635.244116
	n0i[18][3] = 4.00392
	n0i[18][4] = 0.01059
	n0i[18][5] = 0.98763
	n0i[18][6] = 3.06904
	n0i[18][7] = 0
	n0i[18][1] = 21.57882705
	n0i[18][2] = -7766.733078
	n0i[19][3] = 4
	n0i[19][4] = 3.11942
	n0i[19][5] = 1.00243
	n0i[19][6] = 0
	n0i[19][7] = 0
	n0i[19][1] = 21.5830944
	n0i[19][2] = -6069.035869
	n0i[20][3] = 2.5
	n0i[20][4] = 0
	n0i[20][5] = 0
	n0i[20][6] = 0
	n0i[20][7] = 0
	n0i[20][1] = 10.04639507
	n0i[20][2] = -745.375
	n0i[21][3] = 2.5
	n0i[21][4] = 0
	n0i[21][5] = 0
	n0i[21][6] = 0
	n0i[21][7] = 0
	n0i[21][1] = 10.04639507
	n0i[21][2] = -745.375
	th0i[1][4] = 820.659
	th0i[1][5] = 178.41
	th0i[1][6] = 1062.82
	th0i[1][7] = 1090.53
	th0i[2][4] = 662.738
	th0i[2][5] = 680.562
	th0i[2][6] = 1740.06
	th0i[2][7] = 0
	th0i[3][4] = 919.306
	th0i[3][5] = 865.07
	th0i[3][6] = 483.553
	th0i[3][7] = 341.109
	th0i[4][4] = 559.314
	th0i[4][5] = 223.284
	th0i[4][6] = 1031.38
	th0i[4][7] = 1071.29
	th0i[5][4] = 479.856
	th0i[5][5] = 200.893
	th0i[5][6] = 955.312
	th0i[5][7] = 1027.29
	th0i[6][4] = 438.27
	th0i[6][5] = 198.018
	th0i[6][6] = 1905.02
	th0i[6][7] = 893.765
	th0i[7][4] = 468.27
	th0i[7][5] = 183.636
	th0i[7][6] = 1914.1
	th0i[7][7] = 903.185
	th0i[8][4] = 292.503
	th0i[8][5] = 910.237
	th0i[8][6] = 1919.37
	th0i[8][7] = 0
	th0i[9][4] = 178.67
	th0i[9][5] = 840.538
	th0i[9][6] = 1774.25
	th0i[9][7] = 0
	th0i[10][4] = 182.326
	th0i[10][5] = 859.207
	th0i[10][6] = 1826.59
	th0i[10][7] = 0
	th0i[11][4] = 169.789
	th0i[11][5] = 836.195
	th0i[11][6] = 1760.46
	th0i[11][7] = 0
	th0i[12][4] = 158.922
	th0i[12][5] = 815.064
	th0i[12][6] = 1693.07
	th0i[12][7] = 0
	th0i[13][4] = 156.854
	th0i[13][5] = 814.882
	th0i[13][6] = 1693.79
	th0i[13][7] = 0
	th0i[14][4] = 164.947
	th0i[14][5] = 836.264
	th0i[14][6] = 1750.24
	th0i[14][7] = 0
	th0i[15][4] = 228.734
	th0i[15][5] = 326.843
	th0i[15][6] = 1651.71
	th0i[15][7] = 1671.69
	th0i[16][4] = 2235.71
	th0i[16][5] = 1116.69
	th0i[16][6] = 0
	th0i[16][7] = 0
	th0i[17][4] = 1550.45
	th0i[17][5] = 704.525
	th0i[17][6] = 0
	th0i[17][7] = 0
	th0i[18][4] = 268.795
	th0i[18][5] = 1141.41
	th0i[18][6] = 2507.37
	th0i[18][7] = 0
	th0i[19][4] = 1833.63
	th0i[19][5] = 847.181
	th0i[19][6] = 0
	th0i[19][7] = 0
	th0i[20][4] = 0
	th0i[20][5] = 0
	th0i[20][6] = 0
	th0i[20][7] = 0
	th0i[21][4] = 0
	th0i[21][5] = 0
	th0i[21][6] = 0
	th0i[21][7] = 0

	// Precalculations of constants
	for i := 1; i <= MaxFlds; i++ {
		Ki25[i] = math.Pow(Ki[i], 2.5)
		Ei25[i] = math.Pow(Ei[i], 2.5)
	}
	for i := 1; i <= MaxFlds; i++ {
		for j := i; j <= MaxFlds; j++ {
			for n := 1; n <= 18; n++ {
				Bsnij = 1
				if gn[n] == 1 {
					Bsnij = Gij[i][j] * (Gi[i] + Gi[j]) / 2
				}
				if qn[n] == 1 {
					Bsnij = Bsnij * Qi[i] * Qi[j]
				}
				if fn[n] == 1 {
					Bsnij = Bsnij * Fi[i] * Fi[j]
				}
				if sn[n] == 1 {
					Bsnij = Bsnij * Si[i] * Si[j]
				}
				if wn[n] == 1 {
					Bsnij = Bsnij * Wi[i] * Wi[j]
				}
				Bsnij2[i][j][n] = an[n] * math.Pow(Eij[i][j]*math.Sqrt(Ei[i]*Ei[j]), un[n]) * math.Pow(Ki[i]*Ki[j], 1.5) * Bsnij
			}
			Kij5[i][j] = (math.Pow(Kij[i][j], 5) - 1) * Ki25[i] * Ki25[j]
			Uij5[i][j] = (math.Pow(Uij[i][j], 5) - 1) * Ei25[i] * Ei25[j]
			Gij5[i][j] = (Gij[i][j] - 1) * (Gi[i] + Gi[j]) / 2
		}
	}
	// Ideal gas terms
	d0 = 101.325 / RDetail / 298.15
	for i := 1; i <= MaxFlds; i++ {
		n0i[i][3] = n0i[i][3] - 1
		n0i[i][1] = n0i[i][1] - math.Log(d0)
	}
}

func AGA8(T, P float64, _x []float64) (mm, D, Z, dPdD, dPdD2, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa float64) {
	SetupDetail()
	// var _x = []float64{0.965, 0.003, 0.006, 0.018, 0.0045, 0.001, 0.001, 0.0005, 0.0003, 0.0007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	x := _x[:NcDetail]
	x = append([]float64{0.0}, x...)

	// var mm float64 = 0
	MolarMassDetail(x, &mm)

	var ierr int = 0
	var herr string

	// var (
	// 	T float64 = 280
	// 	P float64 = 6000
	// 	D float64 = 1e10
	// )
	DensityDetail(T, P, x, &D, &ierr, &herr)

	// Sub PropertiesDetail(T, D, x, P, Z, dPdD, dPdD2, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa, A)
	var d2PdTD float64
	// var Z, dPdD, dPdD2, d2PdTD, dPdT, U, H, S, Cv, Cp, W, G, JT, Kappa float64
	PropertiesDetail(T, D, x, &P, &Z, &dPdD, &dPdD2, &d2PdTD, &dPdT, &U, &H, &S, &Cv, &Cp, &W, &G, &JT, &Kappa)

	// fmt.Printf("Inputs-----\n")
	// fmt.Printf("Temperature [K]:                    %0.16g\n", T)
	// fmt.Printf("Pressure [kPa]:                     %0.16g\n", P)
	// fmt.Printf("Outputs-----\n")
	// fmt.Printf("Molar mass [g/mol]:                 %0.16g\n", mm)
	// fmt.Printf("Molar density [mol/l]:              %0.16g\n", D)
	// fmt.Printf("Pressure [kPa]:                     %0.16g\n", P)
	// fmt.Printf("Compressibility factor:             %0.16g\n", Z)
	// fmt.Printf("d(P)/d(rho) [kPa/(mol/l)]:          %0.16g\n", dPdD)
	// fmt.Printf("d^2(P)/d(rho)^2 [kPa/(mol/l)^2]:    %0.16g\n", dPdD2)
	// fmt.Printf("d(P)/d(T) [kPa/K]:                  %0.16g\n", dPdT)
	// fmt.Printf("Energy [J/mol]:                     %0.16g\n", U)
	// fmt.Printf("Enthalpy [J/mol]:                   %0.16g\n", H)
	// fmt.Printf("Entropy [J/mol-K]:                  %0.16g\n", S)
	// fmt.Printf("Isochoric heat capacity [J/mol-K]:  %0.16g\n", Cv)
	// fmt.Printf("Isobaric heat capacity [J/mol-K]:   %0.16g\n", Cp)
	// fmt.Printf("Speed of sound [m/s]:               %0.16g\n", W)
	// fmt.Printf("Gibbs energy [J/mol]:               %0.16g\n", G)
	// fmt.Printf("Joule-Thomson coefficient [K/kPa]:  %0.16g\n", JT)
	// fmt.Printf("Isentropic exponent:                %0.16g\n", Kappa)
	return
}

func main() {
	a := app.New()
	w := a.NewWindow("AGA8")
	w.Resize(fyne.NewSize(1280, 0))

	T := binding.NewFloat()
	P := binding.NewFloat()

	CH4 := binding.NewFloat()
	N2 := binding.NewFloat()
	CO2 := binding.NewFloat()
	C2H6 := binding.NewFloat()
	C3H8 := binding.NewFloat()
	H2O := binding.NewFloat()
	H2S := binding.NewFloat()
	H2 := binding.NewFloat()
	CO := binding.NewFloat()
	O2 := binding.NewFloat()
	iC4H10 := binding.NewFloat()

	nC4H10 := binding.NewFloat()
	iC5H12 := binding.NewFloat()
	nC5H12 := binding.NewFloat()
	nC6H14 := binding.NewFloat()
	nC7H16 := binding.NewFloat()
	nC8H18 := binding.NewFloat()
	nC9H20 := binding.NewFloat()
	nC10H22 := binding.NewFloat()
	He := binding.NewFloat()
	Ar := binding.NewFloat()

	Mm := binding.NewFloat()
	D := binding.NewFloat()
	Z := binding.NewFloat()
	DPdD := binding.NewFloat()
	DPdD2 := binding.NewFloat()
	DPdT := binding.NewFloat()
	U := binding.NewFloat()
	H := binding.NewFloat()
	S := binding.NewFloat()
	Cv := binding.NewFloat()
	Cp := binding.NewFloat()
	W := binding.NewFloat()
	G := binding.NewFloat()
	JT := binding.NewFloat()
	Kappa := binding.NewFloat()

	tempreture := widget.NewEntryWithData(binding.FloatToString(T))
	absolutePressure := widget.NewEntryWithData(binding.FloatToString(P))

	methane := widget.NewEntryWithData(binding.FloatToString(CH4))         // 甲烷
	nitrogen := widget.NewEntryWithData(binding.FloatToString(N2))         // 氮气
	carbonDioxide := widget.NewEntryWithData(binding.FloatToString(CO2))   // 二氧化碳
	ethane := widget.NewEntryWithData(binding.FloatToString(C2H6))         // 乙烷
	propane := widget.NewEntryWithData(binding.FloatToString(C3H8))        // 丙烷
	water := widget.NewEntryWithData(binding.FloatToString(H2O))           // 水
	hydrogenSulfide := widget.NewEntryWithData(binding.FloatToString(H2S)) // 硫化氢
	hydrogen := widget.NewEntryWithData(binding.FloatToString(H2))         // 氢气
	carbonMonoxide := widget.NewEntryWithData(binding.FloatToString(CO))   // 一氧化碳
	oxygen := widget.NewEntryWithData(binding.FloatToString(O2))           // 氧气
	isobutane := widget.NewEntryWithData(binding.FloatToString(iC4H10))    // 异丁烷
	nButane := widget.NewEntryWithData(binding.FloatToString(nC4H10))      // 正丁烷
	isopentane := widget.NewEntryWithData(binding.FloatToString(iC5H12))   // 异戊烷
	nPentane := widget.NewEntryWithData(binding.FloatToString(nC5H12))     // 正戊烷
	nHexane := widget.NewEntryWithData(binding.FloatToString(nC6H14))      // 正己烷
	nHeptane := widget.NewEntryWithData(binding.FloatToString(nC7H16))     // 正庚烷
	nOctane := widget.NewEntryWithData(binding.FloatToString(nC8H18))      // 正辛烷
	nNonane := widget.NewEntryWithData(binding.FloatToString(nC9H20))      // 正壬烷
	nDecane := widget.NewEntryWithData(binding.FloatToString(nC10H22))     // 正癸烷
	helium := widget.NewEntryWithData(binding.FloatToString(He))           // 氦气
	argon := widget.NewEntryWithData(binding.FloatToString(Ar))            // 氩气

	MolarMass := widget.NewLabelWithData(binding.FloatToStringWithFormat(Mm, "%0.16g"))
	MolarDensity := widget.NewLabelWithData(binding.FloatToStringWithFormat(D, "%0.16g"))
	Pressure := widget.NewLabelWithData(binding.FloatToStringWithFormat(P, "%0.16g"))
	CompressibilityFactor := widget.NewLabelWithData(binding.FloatToStringWithFormat(Z, "%0.16g"))
	dPdrho := widget.NewLabelWithData(binding.FloatToStringWithFormat(DPdD, "%0.16g"))
	d2Pdrho2 := widget.NewLabelWithData(binding.FloatToStringWithFormat(DPdD2, "%0.16g"))
	dPdT := widget.NewLabelWithData(binding.FloatToStringWithFormat(DPdT, "%0.16g"))
	Energy := widget.NewLabelWithData(binding.FloatToStringWithFormat(U, "%0.16g"))
	Enthalpy := widget.NewLabelWithData(binding.FloatToStringWithFormat(H, "%0.16g"))
	Entropy := widget.NewLabelWithData(binding.FloatToStringWithFormat(S, "%0.16g"))
	IsochoricHeatCapacity := widget.NewLabelWithData(binding.FloatToStringWithFormat(Cv, "%0.16g"))
	IsobaricHeatCapacity := widget.NewLabelWithData(binding.FloatToStringWithFormat(Cp, "%0.16g"))
	SpeedOfSound := widget.NewLabelWithData(binding.FloatToStringWithFormat(W, "%0.16g"))
	GibbsEnergy := widget.NewLabelWithData(binding.FloatToStringWithFormat(G, "%0.16g"))
	JouleThomsonCoefficient := widget.NewLabelWithData(binding.FloatToStringWithFormat(JT, "%0.16g"))
	IsentropicExponent := widget.NewLabelWithData(binding.FloatToStringWithFormat(Kappa, "%0.16g"))

	cleft := container.New(
		layout.NewFormLayout(),
		widget.NewLabel("CH4"), methane,
		widget.NewLabel("N2"), nitrogen,
		widget.NewLabel("CO2"), carbonDioxide,
		widget.NewLabel("C2H6"), ethane,
		widget.NewLabel("C3H8"), propane,
		widget.NewLabel("H2O"), water,
		widget.NewLabel("H2S"), hydrogenSulfide,
		widget.NewLabel("H2"), hydrogen,
		widget.NewLabel("CO"), carbonMonoxide,
		widget.NewLabel("O2"), oxygen,
		widget.NewLabel("i-C4H10"), isobutane,
	)

	cright := container.New(
		layout.NewFormLayout(),
		widget.NewLabel("n-C4H10"), nButane,
		widget.NewLabel("i-C5H12"), isopentane,
		widget.NewLabel("n-C5H12"), nPentane,
		widget.NewLabel("n-C6H14"), nHexane,
		widget.NewLabel("n-C7H16"), nHeptane,
		widget.NewLabel("n-C8H18"), nOctane,
		widget.NewLabel("n-C9H20"), nNonane,
		widget.NewLabel("n-C10H22"), nDecane,
		widget.NewLabel("He"), helium,
		widget.NewLabel("Ar"), argon,
	)

	composition := container.NewGridWithColumns(
		2,
		cleft, cright,
	)

	condition := container.New(
		layout.NewFormLayout(),
		widget.NewLabel("tempreture [K]"), tempreture,
		widget.NewLabel("absolute pressure [kPa]"), absolutePressure,
	)

	reset := widget.NewButton("Reset", func() {
		T.Set(0)
		P.Set(0)
		CH4.Set(0)
		N2.Set(0)
		CO2.Set(0)
		C2H6.Set(0)
		C3H8.Set(0)
		H2O.Set(0)
		H2S.Set(0)
		H2.Set(0)
		CO.Set(0)
		O2.Set(0)
		iC4H10.Set(0)
		nC4H10.Set(0)
		iC5H12.Set(0)
		nC5H12.Set(0)
		nC6H14.Set(0)
		nC7H16.Set(0)
		nC8H18.Set(0)
		nC9H20.Set(0)
		nC10H22.Set(0)
		He.Set(0)
		Ar.Set(0)

		Mm.Set(0)
		D.Set(0)
		Z.Set(0)
		DPdD.Set(0)
		DPdD2.Set(0)
		DPdT.Set(0)
		U.Set(0)
		H.Set(0)
		S.Set(0)
		Cv.Set(0)
		Cp.Set(0)
		W.Set(0)
		G.Set(0)
		JT.Set(0)
		Kappa.Set(0)
	})
	calculate := widget.NewButton("Calculate", func() {
		t, _ := T.Get()
		p, _ := P.Get()

		ch4, _ := CH4.Get()
		n2, _ := N2.Get()
		co2, _ := CO2.Get()
		c2h6, _ := C2H6.Get()
		c3h8, _ := C3H8.Get()
		h2o, _ := H2O.Get()
		h2s, _ := H2S.Get()
		h2, _ := H2.Get()
		co, _ := CO.Get()
		o2, _ := O2.Get()
		ic4h10, _ := iC4H10.Get()
		nc4h10, _ := nC4H10.Get()
		ic5h12, _ := iC5H12.Get()
		nc5h12, _ := nC5H12.Get()
		nc6h14, _ := nC6H14.Get()
		nc7h16, _ := nC7H16.Get()
		nc8h18, _ := nC8H18.Get()
		nc9h20, _ := nC9H20.Get()
		nc10h22, _ := nC10H22.Get()
		he, _ := He.Get()
		ar, _ := Ar.Get()

		mm, d, z, dpdd, dpdd2, dpdt, u, h, s, cv, cp, w, g, jt, kappa := AGA8(t, p, []float64{ch4, n2, co2, c2h6, c3h8, ic4h10, nc4h10, ic5h12, nc5h12, nc6h14, nc7h16, nc8h18, nc9h20, nc10h22, h2, o2, co, h2o, h2s, he, ar})
		Mm.Set(mm)
		D.Set(d)
		Z.Set(z)
		DPdD.Set(dpdd)
		DPdD2.Set(dpdd2)
		DPdT.Set(dpdt)
		U.Set(u)
		H.Set(h)
		S.Set(s)
		Cv.Set(cv)
		Cp.Set(cp)
		W.Set(w)
		G.Set(g)
		JT.Set(jt)
		Kappa.Set(kappa)
	})
	buttonGroup := container.NewGridWithRows(1, calculate, reset)

	title := container.NewCenter(widget.NewLabelWithStyle("Natural Gas Calculation", fyne.TextAlignCenter, fyne.TextStyle{
		Bold:      true,
		Italic:    true,
		Monospace: true,
	}))
	left := container.NewVBox(widget.NewLabel("Composition"), composition, widget.NewLabel("Condition"), condition, layout.NewSpacer(), buttonGroup)
	output := container.New(
		layout.NewFormLayout(),
		widget.NewLabel("Molar mass [g/mol]:"), MolarMass,
		widget.NewLabel("Molar density [mol/l]:"), MolarDensity,
		widget.NewLabel("Pressure [kPa]:"), Pressure,
		widget.NewLabel("Compressibility factor:"), CompressibilityFactor,
		widget.NewLabel("d(P)/d(rho) [kPa/(mol/l)]:"), dPdrho,
		widget.NewLabel("d^2(P)/d(rho)^2 [kPa/(mol/l)^2]:"), d2Pdrho2,
		widget.NewLabel("d(P)/d(T) [kPa/K]:"), dPdT,
		widget.NewLabel("Energy [J/mol]:"), Energy,
		widget.NewLabel("Enthalpy [J/mol]:"), Enthalpy,
		widget.NewLabel("Entropy [J/mol-K]:"), Entropy,
		widget.NewLabel("Isochoric heat capacity [J/mol-K]:"), IsochoricHeatCapacity,
		widget.NewLabel("Isobaric heat capacity [J/mol-K]:"), IsobaricHeatCapacity,
		widget.NewLabel("Speed of sound [m/s]:"), SpeedOfSound,
		widget.NewLabel("Gibbs energy [J/mol]:"), GibbsEnergy,
		widget.NewLabel("Joule-Thomson coefficient [K/kPa]:"), JouleThomsonCoefficient,
		widget.NewLabel("Isentropic exponent:"), IsentropicExponent,
	)

	right := container.NewVBox(output, layout.NewSpacer())
	main := container.NewVBox(title, container.NewGridWithColumns(2, left, right))

	transT := binding.NewFloat()
	transP := binding.NewFloat()
	transM := binding.NewFloat()

	TransT := widget.NewEntryWithData(binding.FloatToString(transT))
	TransP := widget.NewEntryWithData(binding.FloatToString(transP))
	TransM := widget.NewEntryWithData(binding.FloatToString(transM))

	baseInfo := container.New(
		layout.NewFormLayout(),
		widget.NewLabel("temperature [C]"), TransT,
		widget.NewLabel("absolute pressure [Pa]"), TransP,
		widget.NewLabel("Molar mass [g/mol]"), TransM,
	)

	ppm := binding.NewFloat()
	mgm := binding.NewFloat()

	PPM := widget.NewEntryWithData(binding.FloatToString(ppm))
	MGM := widget.NewEntryWithData(binding.FloatToString(mgm))

	ppmArea := container.New(layout.NewFormLayout(), widget.NewLabel("ppm"), PPM)
	mgmArea := container.New(layout.NewFormLayout(), widget.NewLabel("mg/m^3"), MGM)

	leftBtn := widget.NewButton("<", func() {
		mgm.Set(0)
		tPpm, _ := ppm.Get()
		m, _ := transM.Get()
		t, _ := transT.Get()
		p, _ := transP.Get()
		mgm.Set(tPpm * m * 273.15 * p / (22.4 * (273.15 + t) * 101325))
	})
	rightBtn := widget.NewButton(">", func() {
		ppm.Set(0)
		tMgm, _ := mgm.Get()
		m, _ := transM.Get()
		t, _ := transT.Get()
		p, _ := transP.Get()
		ppm.Set((tMgm * 22.4 * (273.15 + t) * 101325) / (m * 273.15 * p))
	})

	transArea := container.NewGridWithColumns(3,
		mgmArea,
		container.NewCenter(container.NewHBox(leftBtn, rightBtn)),
		ppmArea,
	)

	transform := container.NewVBox(
		widget.NewLabel("Base Info"),
		baseInfo,
		widget.NewLabel("Transform Area"),
		transArea,
	)

	helpB64 := `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`

	source := base64.NewDecoder(base64.StdEncoding, strings.NewReader(helpB64))
	help := canvas.NewImageFromReader(source, "help")
	help.FillMode = canvas.ImageFillOriginal

	tabs := container.NewAppTabs(
		container.NewTabItem("Main", container.NewPadded(main)),
		container.NewTabItem("Transform", container.NewPadded(transform)),
		// container.NewTabItem("Help", container.NewPadded(help)),
		container.NewTabItem("About", container.NewCenter(
			container.NewVBox(
				widget.NewLabel("Author: -LAN-"),
				widget.NewLabel("Github: https://github.com/laipz8200/AGA8"),
			),
		)),
	)
	w.SetContent(container.NewPadded(tabs))
	w.ShowAndRun()
}
