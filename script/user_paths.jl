

computername = get(ENV, "COMPUTERNAME", "")

if computername == "SLFW28283"
    shared_folder = "C:/Users/magnusso/OneDrive - Eidg. Forschungsanstalt WSL/Statkraft"
    local_folder = "D:/snowinflow_project/snowinflow_data"
elseif computername == "SLFW29686"
    shared_folder = "C:/Users/opreljen/OneDrive - Eidg. Forschungsanstalt WSL/Statkraft"
    local_folder = "C:/Users/opreljen/LocalDocuments/snowinflow_data"
elseif computername == "SLFW29838"
    shared_folder = "D:/snowinflow_project/snowinflow_shared"
    local_folder = "D:/snowinflow_project/snowinflow_data"
elseif computername == "N51532X"
    shared_folder = "/home/u50116/projects/snow_inflow"
    local_folder = shared_folder
else
    @error "Unknown computer name: $computername"
    shared_folder = ""
    local_folder = ""
end