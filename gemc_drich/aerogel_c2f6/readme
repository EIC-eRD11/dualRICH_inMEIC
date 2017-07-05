This is a standalone RICH simulation using "eic_gemc"

source your set_eic 

./det1_rich_dual_ionside.pl config_det1_rich_dual_ionside.dat  (geometry, hit, bank file produced)

eic_gemc det1_rich_dual_ionside_inter.gcard   (run simulation in graphic mode)

eic_gemc det1_rich_dual_ionside_batch.gcard   (run simulation in batch mode)

more output.txt            (view text output)

evio2root -INPUTF=output.evio  -B=det1_rich_dual_ionside (convert evio output into root tree)

root output.root           (view root output)
