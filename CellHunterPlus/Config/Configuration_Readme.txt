//==================================================================

// parameters_CellHunterPlus.csv: Configuration file for CellHunter+ parameters

//-----------------------------------------------------------------------------------------------------------------

// The settings here are supplied as an example (settings used for CausalXtract publication) 

// and need to be adapted to your experiment.

// polarity="bright"; flag_imm=1; r_tu=14; Rmax_tu=40; DP_tu=70; dist_tu=30;  
// r_sp=4; Rmax_sp=20; DP_sp=10; L_sp=10; r_std=4

//-----------------------------------------------------------------------------------------------------------------

// Parameters description: 

// "polarity" : it must be set to "bright" if bright cells are identified in a dark background.   
// Otherwise, it must be set to "dark".

// "flag_imm" : it must be set to 1 if you want to consider the presence of immune cells,   
// otherwise it must be set to 0. 

// "r_tu": it is the theoretical radius for detecting cancer cells, in pixels.    

// "Rmax_tu": it is the maximum distance for tracking cancer cells,  
// i.e. for linking two presumed instances of the same cancer cell in two different frames  
// to construct the trajectory of that cancer cell, in pixels.  

// "DP_tu": it is the number of frames after which the trajectory of a cancer cell is stopped   
// if the cancer cell is not detected for that specific number of frames.    

// "dist_tu": it imposes to detect only the cancer cells whose centre is less than  
// "dist_tu" pixels away from the centre of the ROI.  

// "r_sp": it is the theoretical radius for detecting immune cells, in pixels.  

// "Rmax_sp": is the maximum distance for tracking immune cells,  
// i.e. for linking two presumed instances of the same immune cell in two different frames  
// to construct the trajectory of that immune cell, in pixels.  

// "DP_sp": it is the number of frames after which the trajectory of an immune cell is stopped  
// if the immune cell is not detected for that specific number of frames.  

// "L_sp": in frames, it filters the immune cells trajectories that are returned as output  
// in the tracking refining process based on their length which must be less than "L_sp".   

// "r_std" in pixels, allows to delete the trajectories of presumed immune cells
// that do not move enough to be considered as such. 
// If r_std is increased, less immune cells are considered. 
// If r_std is decreased, more immune cells are considered.  

//==================================================================

// global_division.csv: Configuration file for tracking correction

//-----------------------------------------------------------------------------------------------------------------

// The settings here are supplied as an example (settings used for CausalXtract publication) 

// and need to be adapted to your experiment.

// If you do not wish to apply the tracking division correction, delete the global_division.csv file.

//-----------------------------------------------------------------------------------------------------------------

// Column description: 

// "IDExp": file name of the roi.

// "global_division": it must be set to 1 if the Main Cancer Cell (MCC) in the roi undergoes division,  
// otherwise it must be set to 0.

// "frame_division": it must be set to the number of frame where the MCC undergoes division,  
// otherwise it must be set to NaN.

//==================================================================

// state_data.csv: Configuration file for "state" features

//-----------------------------------------------------------------------------------------------------------------

// The settings here are supplied as an example (settings used for CausalXtract publication) 

// and need to be adapted to your experiment.

// If you do not have one condition among "CAF_presence", "treatment", "apoptosis", "division", delete the corresponding column.

//-----------------------------------------------------------------------------------------------------------------

// "IDExp": file name of the roi.  

// "ID_frame": number of frame.  

// "CAF_presence": it must be set to 1 if CAFs are present, otherwise it must be set to 0.  

// "treatment": it must be set to 1 if treatment is present, otherwise it must be set to 0.  

// "apoptosis": it indicates if the cell has died during the experiment.  
// It must be set to 0 as long as the MCC is alive, otherwise it must be set to 1 if the MCC has died.  

// "division": it indicates if the cell has divided during the experiment.  
// It must be set to 0 if the MCC has not divided yet, otherwise it must be set to 1 if the MCC has divided.  
