"""
    searchdir(path, key)

Search for files in a directory that contain a specific key string in their filename.
"""
searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))