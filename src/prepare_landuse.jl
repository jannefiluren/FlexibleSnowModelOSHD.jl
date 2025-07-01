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

  return landuse

end
