def import_with_auto_install(packages, scope=locals()):
    if isinstance(packages, str): packages=[packages]
    for package in packages:
        if package.find(":")!=-1:
            package_import_name, package_pip_name = package.split(":")
        else:
            package_import_name, package_pip_name = package, package
        try:
            scope[package_import_name] = __import__(package_import_name)
        except ImportError:
            import subprocess
            subprocess.call(f'pip install {package_pip_name}', shell=True)
            scope[package_import_name] =  __import__(package_import_name)
            
import sys
import os
import shutil
from pathlib import Path

if Path("/home/appuser").exists():
    # essential to avoid cctbx import errors
    target = Path("/home/appuser/venv/share/cctbx")
    if not target.exists():
        target.symlink_to("/home/appuser/.conda/share/cctbx")

    sys.path += ["/home/appuser/venv/lib/python3.9/lib-dynload"]
    os.environ["PATH"] += os.pathsep + "/home/appuser/.conda/bin" 

import streamlit as st
import numpy as np
import mrcfile
import pandas
import re
from shutil import which
from findmysequence_lib.findmysequence import fms_main
import pickle
from bokeh.plotting import ColumnDataSource, figure, output_file, save
from bokeh.models import Label

tmpdir = "map2seq_out"
if not os.path.isdir(tmpdir):
    os.mkdir(tmpdir)

tophits = "--tophits 30000"
hmmer_out = "hmmer_output.txt"

matches_found = "Best matches"
no_matches_found = "No matches found"

def main():
    title = "Identification of Proteins from Density Map"
    st.set_page_config(page_title=title, layout="wide")
    st.title(title)

    #https://discuss.streamlit.io/t/hide-titles-link/19783/4
    st.markdown(""" <style> .css-15zrgzn {display: none} </style> """, unsafe_allow_html=True)

    with st.expander(label="README", expanded=False):
            st.write("This is a Web App to help users to identify proteins that best explain a density map. The user will provide a density map (mrc or ccp4 map format) and a Calpha model in pdb format. The user can also specify a database of protein sequences for the search. Currently, the [findMySequence](https://journals.iucr.org/m/issues/2022/01/00/pw5018/) method is included although we plan to include additional methods in the future.  \nNOTE: the uploaded map/model files are strictly confidential. The developers of this app does not have access to the files")

    col1, col2 = st.columns([1, 3])

    mrc = None
    pdb = None

    with col1:
        #st.subheader("Settings")
        #st.markdown(which("hmmsearch"))

        mrc = None
        # make radio display horizontal
        #st.markdown('<style>div.row-widget.stRadio > div{flex-direction:row;}</style>', unsafe_allow_html=True)
        # input_modes_map = {0:"upload", 1:"emd-xxxxx"}
        input_modes_map = {0:"upload", 1:"url", 2:"emd-xxxxx"}
        help_map = "Only maps in MRC (*\*.mrc*) or CCP4 (*\*.map*) format are supported. Compressed maps (*\*.gz*) will be automatically decompressed"
        input_mode_map = st.radio(label="How to obtain the input map:", options=list(input_modes_map.keys()), format_func=lambda i:input_modes_map[i], index=2, horizontal=True, help=help_map, key="input_mode_map")
        #is_emd = False
        emdb_ids_all, emdb_ids_helical, methods = get_emdb_ids()
        if input_mode_map == 0: # "upload a MRC file":
                label = "Upload a map in MRC or CCP4 format"
                help = None
                fileobj = st.file_uploader(label, type=['mrc', 'map', 'map.gz'], help=help, key="file_upload")
                if fileobj is not None:
                    remove_old_maps()
                    emd_id = extract_emd_id(fileobj.name)
                    is_emd = emd_id is not None
                    with open(os.path.join(tmpdir, fileobj.name), "wb") as f:
                        f.write(fileobj.getbuffer())
                    mrc = tmpdir + "/" + fileobj.name
        elif input_mode_map == 1: # "url":
            url_default = "https://ftp.wwpdb.org/pub/emdb/structures/EMD-10499/map/emd_10499.map.gz"
            help = "An online url (http:// or ftp://) or a local file path (/path/to/your/structure.mrc)"
            url = st.text_input(label="Input the url of a 3D map:", value=url_default, help=help, key="url")
            emd_id = extract_emd_id(url)
            is_emd = emd_id is not None and emd_id
            with st.spinner(f'Downloading {url.strip()}'):
                mrc = get_file_from_url(url.strip())
        elif input_mode_map == 2: # "emdb": randomly selects form emdb_ids_all not emdb_ids_helical anymore
            if not emdb_ids_all:
                st.warning("failed to obtained a list of helical structures in EMDB")
                return
            url = "https://www.ebi.ac.uk/emdb/search/*%20AND%20structure_determination_method:%22helical%22?rows=10&sort=release_date%20desc"
            #st.markdown(f'[All {len(emdb_ids_helical)} helical structures in EMDB]({url})')
            st.markdown(f'[All {len(emdb_ids_all)} structures in EMDB]({url})')
            emd_id_default = "emd-3488"
            do_random_embid = st.checkbox("Choose a random EMDB ID", value=False, key="random_embid")
            if do_random_embid:
                help = "Randomly select another helical structure in EMDB"
                # if max_map_size>0: help += f". {warning_map_size}"
                button_clicked = st.button(label="Change EMDB ID", help=help)
                if button_clicked:
                    import random
                    st.session_state.emd_id = 'emd-' + random.choice(emdb_ids_all)
            else:
                help = None
                # if max_map_size>0: help = warning_map_size
                label = "Input an EMDB ID (emd-xxxxx):"
                emd_id = st.text_input(label=label, value=emd_id_default, key='emd_id', help=help)
                emd_id = emd_id.lower().split("emd-")[-1]
                if emd_id not in emdb_ids_all:
                    import random
                    msg = f"EMD-{emd_id} is not a valid EMDB entry. Please input a valid id (for example, a randomly selected valid id 'emd-{random.choice(emdb_ids_helical)}')"
                    st.warning(msg)
                    return
                #elif emd_id not in emdb_ids_helical:
                #    msg= f"EMD-{emd_id} is in EMDB but annotated as a '{methods[emdb_ids_all.index(emd_id)]}' structure, not a helical structure" 
                #    st.warning(msg)
            if 'emd_id' in st.session_state: emd_id = st.session_state.emd_id
            else: emd_id = emd_id_default
            emd_id = emd_id.lower().split("emd-")[-1]
            with st.spinner(f'Downloading EMD-{emd_id} from {get_emdb_map_url(emd_id)}'):
                filepath = get_emdb_map(emd_id)
                mrc = filepath
            if mrc is None:
                st.warning(f"Failed to download [EMD-{emd_id}](https://www.ebi.ac.uk/emdb/entry/EMD-{emd_id})")
                return

        if mrc is None:
            st.warning(f"Failed to load density map")
            return

        #pdb input
        input_modes_model = {0:"upload", 1:"url", 2:"PDB ID"}
        #input_modes_model = {0:"upload"}
        help_model = "The input PDB model should have all backbone atoms (C-alpha,N,C) of each residue. Sidechain atoms are not required, resiudes can be labeled as any amino acids."
        input_mode_model = st.radio(label="How to obtain the input PDB file:", options=list(input_modes_model.keys()), format_func=lambda i:input_modes_model[i], index=2, horizontal=True, help=help_model, key="input_mode_model")
        # pdb_ids_all = get_pdb_ids()
        pdb = None

        if input_mode_model == 0: # "upload a PDB file":
            label = "Upload a PDB file"
            fileobj = st.file_uploader(label, type=['pdb'], help=None, key="file_upload")
            if fileobj is not None:
                remove_old_pdbs()
                with open(os.path.join(tmpdir, fileobj.name), "wb") as f:
                    f.write(fileobj.getbuffer())
                pdb = tmpdir + "/" + fileobj.name
                # pdb = get_model_from_uploaded_file(fileobj)
        
        elif input_mode_model == 1: # "url":
            help = "An online url (http:// or ftp://) or a local file path (/path/to/your/model.pdb)"
            url = st.text_input(label="Input the url of a PDB model:", help=help, key="url")
            if url:
                with st.spinner(f'Downloading {url.strip()}'):
                    remove_old_pdbs()
                    pdb = get_file_from_url(url.strip())
        
        elif input_mode_model == 2: # "PDB ID":
            help = None
            # if max_map_size>0: help = warning_map_size
            label = "Input an PDB ID (for example: 4hhb):"
            pdb_id_default = "5me2"
            pdb_id = st.text_input(label=label, key='pdb_id', value=pdb_id_default, help=help)
            pdb_id = pdb_id.lower()
            if pdb_id:
                pdb_url=get_pdb_url(pdb_id)
                with st.spinner(f'Downloading {pdb_id}.pdb from {pdb_url}'):
                    pdb = get_file_from_url(pdb_url)
        
        if pdb is None:
            st.warning(f"Failed to load the PDB model")
            return

        #db input
        input_modes_db = {0:"upload", 1:"url", 2:"human proteins"}
        #input_modes_model = {0:"upload"}
        help_db = "The input sequence library in .fa.gz"
        input_mode_db = st.radio(label="Which sequence database to use:", options=list(input_modes_db.keys()), format_func=lambda i:input_modes_db[i], index=2, horizontal=True, help=help_db, key="input_mode_db")
        # pdb_ids_all = get_pdb_ids()
        
        db = None

        if input_mode_db == 0: # "upload":
            label = "Upload a compressed fasta file (.fa.gz)"
            fileobj = st.file_uploader(label, type=['fa.gz'], help=None, key="file_upload")
            if fileobj is not None:
                #remove_old_db()
                with open(os.path.join(tmpdir, fileobj.name), "wb") as f:
                    f.write(fileobj.getbuffer())
                db = tmpdir + "/" + fileobj.name
                # pdb = get_model_from_uploaded_file(fileobj)
        else:
            if input_mode_db == 1: # "url":
                help = "An online url (http:// or ftp://) or a local file path (/path/to/your/database.fa.gz)"
                url = st.text_input(label="Input the url of a sequence database (.fa.gz):", help=help, key="url")
            elif input_mode_db == 2: # "use default":
                url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
                st.markdown(f"Using the [human protein sequences]({url})")
            with st.spinner(f'Downloading {url.strip()}'):
                #remove_old_db()
                db = get_file_from_url(url.strip())
        
        if db is None:
            st.warning(f"Failed to load the protein sequence database")
            return
        
        direction_options = {0:"original", 1:"reversed"}
        help_direction=None
        direction_option = st.radio(label="Protein sequence direction:", options=list(direction_options.keys()), format_func=lambda i:direction_options[i], index=0, horizontal=True, help=help_direction, key="direction_option")
        
        handedness_options = {0:"original", 1:"flipped"}
        help_handedness=None
        handedness_option = st.radio(label="Map handedness:", options=list(handedness_options.keys()), format_func=lambda i:handedness_options[i], index=0, horizontal=True, help=help_handedness, key="handedness_option")
        
        if handedness_option in [1, 2]: # flipped
            flip_map_model(mrc, pdb)    

    with col2:
        with st.spinner("Processing..."):
            #remove_old_graph_log()
            seqin = None
            modelout = None
            res = map2seq_run(mrc, pdb, seqin, modelout, direction_option, handedness_option, db, outdir = tmpdir)
            if res is None:
                st.error(f"Failed")
                return

            xs, ys = res
                
        #https://docs.bokeh.org/en/latest/docs/user_guide/tools.html
        
        source = ColumnDataSource(data=dict(x=range(len(xs)),y=ys,ID=xs))
        top_source = ColumnDataSource(data=dict(x=[0],y=[ys[0]],ID=[xs[0]]))
        label = Label(x=0, y=ys[0], text='Best Match', x_offset=10, y_offset=-5, render_mode='canvas')

        TOOLTIPS = [('index','$index'),('ID','@ID'),('E-val','@y')]

        p = figure(tooltips=TOOLTIPS, y_axis_type='log', title='Ranked Sequences')
        p.circle('x','y',source=source)
        p.circle('x','y',source=top_source, size=10,line_color='red',fill_color='red')
        p.yaxis.axis_label = 'E-values'
        p.xaxis.axis_label = 'Rank Order'
        p.y_range.flipped = True
        p.add_layout(label)
        
        st.bokeh_chart(p, use_container_width=True)
        
        ## Prepare for the second run
        #import pyfastx
        #fa = pyfastx.Fasta("./tempDir/human.fa.gz")
        #seqin = tmpdir+"/tmp.fasta"
        #modelout = tmpdir+"/model_out.pdb"
        #with open(tmpdir+"/tmp.fasta","w") as tmp:
        #    tmp.write(">"+xs[0]+"\n")
        #    tmp.write(fa[xs[0]].seq)
        #
        #map2seq_run(mrc, pdb, seqin, modelout, direction_option, handedness_option, db, outdir = tmpdir)
        #
        #st.write("Alignment with "+xs[0]+":")
        #
        #with open(tmpdir+"/seq_align_output.txt","r") as tmp:
        #    for line in tmp.readlines():
        #        if line[0:7]!="WARNING":
        #            st.write(line)
        #            
        #with open(modelout,"r") as tmp:
        #    out_texts="".join(tmp.readlines())
        #    st.download_button("Download output model", data=out_texts, file_name="model_out.pdb")
        
    with col2:
        df = pandas.DataFrame({"E-val (log10)":np.log10(ys).T, "Protein":np.array(xs).T})
        df.index += 1
              
        n = 10
        st.subheader(f"Top {n} matches:")
        df_top = df.iloc[:n, :].copy()
        def link_to_uniprot(s):
            pid = s.split('|')[1]
            url = f"https://www.uniprot.org/uniprotkb/{pid}"
            return f'<a target="_blank" href="{url}">{pid}</a>'            
        df_top.loc[:, "Link"] = df_top.loc[:, "Protein"].apply(link_to_uniprot)
        df_top.loc[:, "Protein"] = df_top.loc[:, "Protein"].str.split("|", expand=True).iloc[:, -1]
        df_top.reset_index(inplace=True)
        df_top = df_top.rename(columns = {'index':'Rank'})
        st.write(df_top.to_html(escape=False, index=False, justify="left"), unsafe_allow_html=True)

        def link_to_uniprot_2(s):
            pid = s.split('|')[1]
            url = f"https://www.uniprot.org/uniprotkb/{pid}"
            return url
        df.loc[:, "Uniprot ID"] = df.loc[:, "Protein"].str.split("|", expand=True).iloc[:, 1]
        df.loc[:, "URL"] = df.loc[:, "Protein"].apply(link_to_uniprot_2)
        df.loc[:, "Protein"] = df.loc[:, "Protein"].str.split("|", expand=True).iloc[:, -1]
        st.download_button(
            label=f"Download the scores for all {len(df)} proteins",
            data=df.to_csv().encode('utf-8'),
            file_name='map2seq_results.csv',
            mime='text/csv'
        )

    #remove_old_pdbs()
    #remove_old_maps()
    #remove_old_graph_log()
    mrc=None
    pdb=None

    
def remove_old_graph_log():
    dir = os.listdir(tmpdir)
    for item in dir:
        if item.endswith(".png") or item.endswith(".txt"):
            os.remove(os.path.join(tmpdir, item))

#-------------------------------Map Functions-------------------------------
def remove_old_maps():

    dir = os.listdir(tmpdir)
    for item in dir:
        if item.endswith(".mrc") or item.endswith(".map") or item.endswith(".map.gz"):
            os.remove(os.path.join(tmpdir, item))

def get_direct_url(url):
    import re
    if url.startswith("https://drive.google.com/file/d/"):
        hash = url.split("/")[5]
        return f"https://drive.google.com/uc?export=download&id={hash}"
    elif url.startswith("https://app.box.com/s/"):
        hash = url.split("/")[-1]
        return f"https://app.box.com/shared/static/{hash}"
    elif url.startswith("https://www.dropbox.com"):
        if url.find("dl=1")!=-1: return url
        elif url.find("dl=0")!=-1: return url.replace("dl=0", "dl=1")
        else: return url+"?dl=1"
    elif url.find("sharepoint.com")!=-1 and url.find("guestaccess.aspx")!=-1:
        return url.replace("guestaccess.aspx", "download.aspx")
    elif url.startswith("https://1drv.ms"):
        import base64
        data_bytes64 = base64.b64encode(bytes(url, 'utf-8'))
        data_bytes64_String = data_bytes64.decode('utf-8').replace('/','_').replace('+','-').rstrip("=")
        return f"https://api.onedrive.com/v1.0/shares/u!{data_bytes64_String}/root/content"
    else:
        return url

#"https://ftp.wwpdb.org/" -> threw a IsADirectoryError: This app has encountered an error. ...
#Since I delete maps from the "tempDir" folder 
@st.experimental_memo(persist='disk', max_entries=1, ttl=60*60, show_spinner=False, suppress_st_warning=True)
def get_file_from_url(url):
    url_final = get_direct_url(url)    # convert cloud drive indirect url to direct url
    ds = np.DataSource(None)
    if not ds.exists(url_final):
        st.error(f"ERROR: {url} could not be downloaded. If this url points to a cloud drive file, make sure the link is a direct download link instead of a link for preview")
        st.stop()
    filename_final = None
    localfile_path  = None
    with ds.open(url) as fp:
        # data = get_3d_map_from_file(fp.name)
        if fp.name.endswith(".gz"):
            filename_final = fp.name[:-3]
            import gzip, shutil
            with gzip.open(fp.name, 'r') as f_in, open(filename_final, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            filename_final = fp.name
  
        import shutil
        shutil.copy2(filename_final, tmpdir)
        localfile_path = Path(tmpdir) / Path(filename_final).name
    return str(localfile_path.resolve())

def extract_emd_id(text):
    import re
    pattern = '.*emd_([0-9]*)\.map.*'

    match = re.search(pattern, text, re.IGNORECASE)
    if match:
        emd_id = match.group(1)
    else:
        emd_id = None
    return emd_id

@st.experimental_memo(persist='disk', max_entries=1, ttl=60*60*24, show_spinner=False, suppress_st_warning=True)
def get_emdb_ids():
    try:
        import_with_auto_install(["pandas"])
        import pandas as pd
        entries_all = pd.read_csv('https://www.ebi.ac.uk/emdb/api/search/current_status:"REL"?wt=csv&download=true&fl=emdb_id,structure_determination_method,resolution')
        methods = list(entries_all["structure_determination_method"])
        entries_helical = entries_all[entries_all["structure_determination_method"]=="helical"]
        emdb_ids_all     = list(entries_all.iloc[:,0].str.split('-', expand=True).iloc[:, 1].values)
        emdb_ids_helical = list(entries_helical.iloc[:,0].str.split('-', expand=True).iloc[:, 1].values)
    except:
        emdb_ids_all = []
        emdb_ids_helical = []
        methods = {}
    return emdb_ids_all, emdb_ids_helical, methods

def get_emdb_map_url(emdid):
    emdid_number = emdid.lower().split("emd-")[-1]
    server = "https://ftp.wwpdb.org/pub"    # Rutgers University, USA
    #server = "https://ftp.ebi.ac.uk/pub/databases" # European Bioinformatics Institute, England
    #server = "http://ftp.pdbj.org/pub" # Osaka University, Japan
    url = f"{server}/emdb/structures/EMD-{emdid_number}/map/emd_{emdid_number}.map.gz"
    return url

def get_pdb_url(protid):
	server = "https://files.rcsb.org/download"
	return f"{server}/{protid}.pdb.gz"
	
@st.experimental_memo(persist='disk', max_entries=1, show_spinner=False, suppress_st_warning=True)
def get_emdb_map(emdid):
    url = get_emdb_map_url(emdid)
    mapfile = get_file_from_url(url)
    # If no this line the returning path will be wrong. WHY???
    return mapfile

#-------------------------------End Map Functions-------------------------------

#-------------------------------Model Functions-------------------------------
def remove_old_pdbs():
    dir = os.listdir(tmpdir)
    for item in dir:
        if item.endswith(".pdb"):
            os.remove(os.path.join(tmpdir, item))

@st.experimental_memo(persist='disk', max_entries=1, ttl=60*60*24*7, show_spinner=False, suppress_st_warning=True)
def get_pdb_ids():
    try:
        url = "ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx"
        ds = np.DataSource(None)
        with ds.open(url) as fp:
            pdb_ids = [line[:4] for line in fp.readlines()[2:] if len(line) > 4]
    except:
        pdb_ids = None
    return pdb_ids
#-------------------------------End Model Functions-------------------------------

def flip_map_model(map_name,pdb_name):
    mrc_data = mrcfile.open(map_name, 'r+')
    v_size=mrc_data.voxel_size
    nx=mrc_data.header['nx']
    ny=mrc_data.header['ny']
    nz=mrc_data.header['nz']
    apix=v_size['z']

    data=mrc_data.data
    new_data=np.empty((nx,ny,nz),dtype=np.float32)
    for i in range(nz):
        new_data[i,:,:]=data[nz-1-i,:,:]
        
    mrc_data.set_data(new_data)
    mrc_data.close()
    
    with open(pdb_name,'r') as f:
        lines=f.readlines()
        orig_atom_lines = []
        for line in lines:
            if line[0:4] == "ATOM":
                orig_atom_lines.append(line)
            elif line[0:3] == "TER":
               orig_atom_lines.append(line)
    with open(pdb_name, 'w') as o:
        for line in orig_atom_lines:
            if line[0:3] != "TER":
                coord_vec = np.array(line[30:54].split()).astype(np.float32)
                new_z=(nz-1)*apix-coord_vec[2]
                line = list(line)
                line[30:54] = list(" " + f'{coord_vec[0]:>7.3f}' + " " + f'{coord_vec[1]:>7.3f}' + " " + f'{new_z:>7.3f}')
                line = "".join(line)
            o.write(line)

#------------------------------- Main Functions-------------------------------
@st.experimental_memo(persist='disk', max_entries=1, ttl=60*60, show_spinner=False, suppress_st_warning=True)
def map2seq_run(map, pdb, seqin, modelout, rev, flip, db, outdir = "tempDir/"):

    map = os.path.abspath(map)
    pdb = os.path.abspath(pdb)
    db = os.path.abspath(db)
    outdir = os.path.abspath(outdir)

    if outdir[-1] != "/":
        outdir += "/"

    basename = f"{Path(map).stem}_{Path(pdb).stem}"

    fms_main.fms_run(mapin=map, modelin=pdb, seqin=seqin, modelout=modelout, db=db, tmpdir=outdir, outdir=outdir, rev=rev, flip=flip, tophits=np.iinfo(np.uint32).max)
    
    #graph fms output
    num = parse_file(f"{outdir}{basename}.png", f"{outdir}{hmmer_out}")
    if num == -1:
        return None # failed

    with open(os.path.join(outdir, f'{basename}.png_x.pkl'),'rb') as inf:
        xs = pickle.load(inf)
    with open(os.path.join(outdir, f'{basename}.png_y.pkl'),'rb') as inf:
        ys = pickle.load(inf)

    return (xs, ys)
        
def make_graph(ids, e_vals, outputFile):
    
    #https://docs.bokeh.org/en/latest/docs/user_guide/tools.html

    output_file('{}.html'.format(outputFile))

    source = ColumnDataSource(data=dict(x=range(len(ids)),y=e_vals,ID=ids))
    top_source = ColumnDataSource(data=dict(x=[0],y=[e_vals[0]],ID=[ids[0]]))
    label = Label(x=0, y=e_vals[0], text='Best Match', x_offset=10, y_offset=-5, render_mode='canvas')
  
    TOOLTIPS = [('index','$index'),('ID','@ID'),('E-val','@y')]
   
    p = figure(width=400,height=400,tooltips=TOOLTIPS,y_axis_type='log', title='Ranked Sequences')
    p.circle('x','y',source=source)
    p.circle('x','y',source=top_source, size=10,line_color='red',fill_color='red')
    p.yaxis.axis_label = 'E-values'
    p.xaxis.axis_label = 'Rank Order'
    p.y_range.flipped = True
    p.add_layout(label)
    
    #save(p)
    
    with open('{}_x.pkl'.format(outputFile),'wb') as o:
        pickle.dump(ids,o,pickle.HIGHEST_PROTOCOL)
    with open('{}_y.pkl'.format(outputFile),'wb') as o:
        pickle.dump(e_vals,o,pickle.HIGHEST_PROTOCOL)


def parse_file(outputFile, filepath):
    ids = []
    e_vals = []
    with open(filepath, 'r') as file:
        firstline = file.readline().rstrip()
        if firstline.find(matches_found) == -1:
            print(no_matches_found)
            return -1
        for line in file:
            line = line.rstrip()
            list = line.split('|')
            list[0:3] = ["|".join(list[0:3])]
            #list = line.split(' ')
            list[0] = list[0].strip()
            list[1] = list[1].removeprefix('E-value=')
            list[1] = float(list[1])
            ids.append(list[0])
            e_vals.append(list[1])
                        
        make_graph(ids, e_vals, outputFile)
    return 1
#------------------------------- End Main Functions-------------------------------

if __name__ == "__main__":
    main()
