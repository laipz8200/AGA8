env GOOS="windows" GOARCH="amd64" CGO_ENABLED="1" CC="x86_64-w64-mingw32-gcc" $GOPATH/bin/fyne package -os windows -icon asserts/icon.png
$GOPATH/bin/fyne package -os darwin -icon asserts/icon.png
