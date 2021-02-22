#!/usr/bin/env bash

inkscape "logo.svg" --export-filename="logo-128x128.png" -w 128 -h 128
inkscape "logo.svg" --export-filename="logo-600x600.png" -w 600 -h 600

inkscape "logo-title.svg" --export-filename="logo-title-200.png" -h 200
