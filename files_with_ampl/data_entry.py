def data_entry(cohesion, friction_ang , dilantancy , max_nl , min_nl,max_nl_1,max_nl_2,fname, dxf_file  ): 
    data = dict()                           #define dictionary to collect all data
    
    data["friction_angle"] = friction_ang   #friction angle
    data["cohesion"] = cohesion             #cohesion
    data["dilantancy"] =  dilantancy        # dilatancy angle
    data["max_normalload"] = max_nl         # Maximum allowable normal load 
    data["min_normalload"] = min_nl         # minimum allowable normal load 
    data["notepad_files"]  = fname          # text file names (the file in the first position is no longer used in this code so no need to change it while the second one contain the loading data)
    data["autocad_file"]   = dxf_file       #DXF file location
    
    left_limit = [67.37,16.51]              #left bottom limit of the region where edges have different cohesions and friction angles
    right_limit = [69.29,18.325]            #right top limit of the region where edges have different cohesions and friction angles
    
    
    cohesion_reg = cohesion                 #cohesion on edges lying in the region defined above
    friction_reg = friction_ang             #Friction angle on edges lying in the region defined above
    
    data["max_nl_2"] = max_nl_2
    data["max_nl_1"] = max_nl_1
    
    data["cohesion_region"] = cohesion_reg
    data["friction_ang_region"] = friction_reg
    
    fname = ["mechanic_data_ex5blocks.txt","mechanic_data_ex5blocks.txt"] 
    dxf_file = "H:\Python-projects\ex5blocks.dxf" 
    
    return data