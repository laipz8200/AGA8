package main

import (
	"log"

	"github.com/laipz8200/AGA8/pkg/gui"
)

func main() {
	app := gui.CreateApp()
	app.Run()
	log.Println("App Start")
}
