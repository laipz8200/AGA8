package main

import (
	"log"

	"github.com/laipz8200/AGA8/pkg/gui"
)

func main() {
	app := gui.CreateApp()
	log.Println("App Start")
	app.Run()
}
