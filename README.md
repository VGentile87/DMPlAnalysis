# DMPlAnalysis
Analysis tool for plasmon
#
# INSTRUCTIONS (V. GENTILE)
#
# Library compilation
source compile.sh
#
# ENVIRONMENTAL VARIABLES SETTINGS (bashrc)
#
cd /path/to/DMPlAnalysis
gedit setup.sh
export DMPLA="/path/to/DMPlAnalysis"
cd
gedit .bashrc
source /path/to/DMPlAnalysis/setup.sh
# Usage
#
cd work_folder (with dm_tracks.dm.root file or dm_tracks_cl.dm.root)
cpstart
source start_root6.sh (follow the instructions)
root -l debug8_test_grain.root (enjoy the result)
#
# Additional info
#
# setting.mac (optical microscope and cuts settings)
# log_run.txt (info about scanned area and anlyzed area)
