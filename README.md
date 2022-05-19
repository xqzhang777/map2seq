# Project Title

## Description
A web app based on: 
Streamlit (https://streamlit.io/)
findMySequence(https://gitlab.com/gchojnowski/findmysequence)
  findMySequence: a neural-network-based approach for identification of unknown proteins in X-ray crystallography and cryo-EM
  Grzegorz Chojnowski, Adam J. Simpkin, Diego A. Leonardo, Wolfram Seifert-Davila, Dan E. Vivas-Ruiz, Ronan M. Keegan, Daniel J. Rigden
  IUCrJ 9.1 (2022)

## Instructions on how to run
First run
```
pipenv shell
```

Run conda deactivate
```
conda deactivate
```

Then start findmysequence env
```
conda activate findmysequence
```

Then start streamlit app
```
streamlit run driver.py
```

## Organization
### driver.py
Contains streamlit app

### main.py
Used for running findmysequence and graphing the output

### graph.py
Graphs the hmmer output

### fms.py
Runs findmysequence with correct filepaths

### tempDir
Contains the database, all map files and pdb files that are uploaded will be stored here and are not deleted (explained in Issues) so delete them whenever after you stop the app

## Changes made in some functions
- Get random emdb id now gets id from emdb_ids_all and not emdb_ids_helical
- Does not check or limit map size anymore
- get_3d_map_from_url() does not return numpy array anymore it instead saves mrc file to "tempDir" and returns a file path to the file it saved

## Issues
### get_3d_map_from_url() (and get_model_from_url())
I modified it to return take file from the temporary directory and copy it into "tempDir" and return the file path to that file in "tempDir".
This function has @st.experimental_singleton() so the return value is stored in cache, and it won't redownload files unless cache is cleared. So that means the map files are also being stored and not deleted. I think this might cause some memory limit issues since map files is not being deleted unless the entire streamlit app is restarted or cache is cleared. For now the map files will be stored in Web_App/tempDir so they must be manually deleted. There is a function remove_old_map() which will delete them so it should be added to main right after when (if?) you decide to manually clear the cache. Haven't implemneted yet, but there will be a similar issue with get_model_from_url().
