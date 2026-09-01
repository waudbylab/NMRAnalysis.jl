"""
prepare_outputfolder(outputfolder)

    This function checks if the outputfolder exists. If it doesn't, it creates it.
    If it does exist, it deletes the folder and its contents, then recreates it.
    Please note that this function will delete all files and subdirectories in the
    specified outputfolder, so use it with caution.
"""
function prepare_outputfolder(outputfolder)
    if !isdir(outputfolder)
        mkdir(outputfolder)
    else
        rm(outputfolder; recursive=true, force=true)
        mkdir(outputfolder)
    end
end