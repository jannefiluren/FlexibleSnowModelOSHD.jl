function prepare_landuse_stations(station::String = "")

  landuse = matread("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_STAT.mat")

  if isempty(station)
    landuse["is_domain"] = ones(Bool, size(landuse["dem"]["data"]))
  else
    landuse["is_domain"] = landuse["acro"] .== station
  end

  landuse["dhxydy"] = ones(size(landuse["dem"]["data"]))
  landuse["sddem"] = ones(size(landuse["dem"]["data"]))
  landuse["Ld"] = ones(size(landuse["dem"]["data"]))      # TODO this is set to zero in the wrapper...

  landuse["slopemu"] = sqrt.((landuse["dhxydy"]./2));
  landuse["xi"] = (sqrt(2)*landuse["sddem"])./landuse["slopemu"];

  return landuse

end


function prepare_landuse_grid()

  landuse = matread("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat")

  landuse["is_domain"] = ones(Bool, size(landuse["dem"]["data"]))

  landuse["dhdxdy"] = landuse["dhdxdy"]["data"]
  landuse["sddem"] = landuse["sd"]["data"]
  landuse["Ld"] = 250*ones(size(landuse["dem"]["data"]))

  landuse["slopemu"] = sqrt.((landuse["dhdxdy"]./2));
  landuse["xi"] = (sqrt(2)*landuse["sddem"])./landuse["slopemu"];

  landuse["x"] = ones(size(landuse["dem"]["data"]))   # TODO hack for adding x
  landuse["y"] = ones(size(landuse["dem"]["data"]))   # TODO hack for adding y

  return landuse

end