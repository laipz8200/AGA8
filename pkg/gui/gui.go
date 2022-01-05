package gui

import (
	"encoding/base64"
	"strings"

	"fyne.io/fyne/v2"
	"fyne.io/fyne/v2/app"
	"fyne.io/fyne/v2/canvas"
	"fyne.io/fyne/v2/container"
	"fyne.io/fyne/v2/data/binding"
	"fyne.io/fyne/v2/layout"
	"fyne.io/fyne/v2/widget"
	"github.com/laipz8200/AGA8/pkg/aga8"
)

func CreateApp() fyne.App {
	app := app.New()
	window := app.NewWindow("AGA8")
	window.Resize(fyne.NewSize(1280, 0))

	title := container.NewCenter(widget.NewLabelWithStyle("Natural Gas Calculation", fyne.TextAlignCenter, fyne.TextStyle{
		Bold:      true,
		Italic:    true,
		Monospace: true,
	}))

	main := MainPage()
	transform := TransformPage()
	// help := ChineseHelpPage()

	tabs := container.NewAppTabs(
		container.NewTabItem("Home", container.NewPadded(main)),
		container.NewTabItem("Transform", container.NewPadded(transform)),
		// container.NewTabItem("Help", container.NewPadded(help)),
		container.NewTabItem("About", container.NewCenter(
			container.NewVBox(
				widget.NewLabel("Author: -LAN-"),
				widget.NewLabel("Github: https://github.com/laipz8200/AGA8"),
			),
		)),
	)
	window.SetContent(
		container.NewVBox(
			title,
			container.NewPadded(tabs),
		),
	)
	window.Show()
	return app
}

func MainPage() *fyne.Container {
	// Conditons
	T := binding.NewFloat()
	P := binding.NewFloat()

	// Compositions
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

	// Results
	density := binding.NewFloat()
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

	// Entries
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

	// Labels
	Density := widget.NewLabelWithData(binding.FloatToStringWithFormat(density, "%0.16g"))
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
		widget.NewLabel("CH4"), methane, widget.NewLabel("N2"), nitrogen, widget.NewLabel("CO2"), carbonDioxide, widget.NewLabel("C2H6"), ethane,
		widget.NewLabel("C3H8"), propane, widget.NewLabel("H2O"), water, widget.NewLabel("H2S"), hydrogenSulfide, widget.NewLabel("H2"), hydrogen,
		widget.NewLabel("CO"), carbonMonoxide, widget.NewLabel("O2"), oxygen, widget.NewLabel("i-C4H10"), isobutane,
	)
	cright := container.New(
		layout.NewFormLayout(),
		widget.NewLabel("n-C4H10"), nButane, widget.NewLabel("i-C5H12"), isopentane, widget.NewLabel("n-C5H12"), nPentane, widget.NewLabel("n-C6H14"), nHexane,
		widget.NewLabel("n-C7H16"), nHeptane, widget.NewLabel("n-C8H18"), nOctane, widget.NewLabel("n-C9H20"), nNonane, widget.NewLabel("n-C10H22"), nDecane,
		widget.NewLabel("He"), helium, widget.NewLabel("Ar"), argon,
	)

	composition := container.NewGridWithColumns(2, cleft, cright)

	condition := container.New(
		layout.NewFormLayout(),
		widget.NewLabel("tempreture [K]"), tempreture, widget.NewLabel("absolute pressure [kPa]"), absolutePressure,
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

		density.Set(0)
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

		mm, d, z, dpdd, dpdd2, dpdt, u, h, s, cv, cp, w, g, jt, kappa := aga8.Calc(t, p, []float64{ch4, n2, co2, c2h6, c3h8, ic4h10, nc4h10, ic5h12, nc5h12, nc6h14, nc7h16, nc8h18, nc9h20, nc10h22, h2, o2, co, h2o, h2s, he, ar})
		density.Set(mm * d)
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
	left := container.NewVBox(widget.NewLabel("Composition"), composition, widget.NewLabel("Condition"), condition, layout.NewSpacer(), buttonGroup)
	output := container.New(
		layout.NewFormLayout(),
		widget.NewLabel("Density [kg/m^3]:"), Density, widget.NewLabel("Molar mass [g/mol]:"), MolarMass, widget.NewLabel("Molar density [mol/l]:"), MolarDensity,
		widget.NewLabel("Pressure [kPa]:"), Pressure, widget.NewLabel("Compressibility factor (Z):"), CompressibilityFactor, widget.NewLabel("d(P)/d(rho) [kPa/(mol/l)]:"), dPdrho,
		widget.NewLabel("d^2(P)/d(rho)^2 [kPa/(mol/l)^2]:"), d2Pdrho2, widget.NewLabel("d(P)/d(T) [kPa/K]:"), dPdT, widget.NewLabel("Energy [J/mol]:"), Energy,
		widget.NewLabel("Enthalpy [J/mol]:"), Enthalpy, widget.NewLabel("Entropy [J/mol-K]:"), Entropy, widget.NewLabel("Isochoric heat capacity [J/mol-K]:"), IsochoricHeatCapacity,
		widget.NewLabel("Isobaric heat capacity [J/mol-K]:"), IsobaricHeatCapacity, widget.NewLabel("Speed of sound [m/s]:"), SpeedOfSound,
		widget.NewLabel("Gibbs energy [J/mol]:"), GibbsEnergy, widget.NewLabel("Joule-Thomson coefficient [K/kPa]:"), JouleThomsonCoefficient,
		widget.NewLabel("Isentropic exponent:"), IsentropicExponent,
	)

	right := container.NewPadded(container.NewVBox(output, layout.NewSpacer()))
	return container.NewGridWithColumns(2, left, right)
}

func TransformPage() *fyne.Container {
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

	return container.NewVBox(widget.NewLabel("Base Info"), baseInfo, widget.NewLabel("Transform Area"), transArea)
}

func ChineseHelpPage() *fyne.Container {
	source := base64.NewDecoder(base64.StdEncoding, strings.NewReader(helpB64))
	help := canvas.NewImageFromReader(source, "help")
	help.FillMode = canvas.ImageFillOriginal
	return container.NewPadded(help)
}
