#!/usr/bin/env bash

# launch the pi-control for a buch of model at once

## hadgem2-es
./1_main_analysis.r h08 hadgem2-es picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r pcr-globwb hadgem2-es picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r watergap2 hadgem2-es picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r lpjml hadgem2-es picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r matsiro hadgem2-es picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r cwatm hadgem2-es picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r CLM45 hadgem2-es picontrol picontrol 1661 1661 1661 1661 yes &

## ipsl
./1_main_analysis.r h08 ipsl-cm5a-lr picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r pcr-globwb ipsl-cm5a-lr picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r watergap2 ipsl-cm5a-lr picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r lpjml ipsl-cm5a-lr picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r matsiro ipsl-cm5a-lr picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r cwatm ipsl-cm5a-lr picontrol picontrol 1661 1661 1661 1661 yes &

## miroc5
./1_main_analysis.r h08 miroc5 picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r pcr-globwb miroc5 picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r watergap2 miroc5 picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r lpjml miroc5 picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r matsiro miroc5 picontrol picontrol 1661 1661 1661 1661 yes &
./1_main_analysis.r cwatm miroc5 picontrol picontrol 1661 1661 1661 1661 yes &
