import sys, os, pickle
from pathlib import Path
import streamlit as st

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
            #if Path("/home/appuser").exists():
            #    subprocess.call(f'/home/appuser/.conda/bin/pip install {package_pip_name}', shell=True)
            #else:
            subprocess.call(f'pip install {package_pip_name}', shell=True)
            scope[package_import_name] =  __import__(package_import_name)


tmpdir = "map2seq_out"
if not os.path.isdir(tmpdir):
    os.mkdir(tmpdir)

#st.info(os.environ["LD_LIBRARY_PATH"])
from shutil import which
#st.info(which("python"))


try:
    import cctbx
    #raise ImportError
except ImportError:
    import zstandard
    import tempfile
    import tarfile
    import numpy as np
   
    #st.info("downloading cctbx-base")
    ds = np.DataSource(tmpdir)
    url_final = "https://drive.google.com/uc?export=download&id=1pWpLoyUOXqTbktqOJ24X5bdXa8u7lb0Y"
    if not ds.exists(url_final):
        print("download error")
    with ds.open(url_final) as fp:
        filename_final = fp.name
    filepath_final = Path(filename_final).resolve()
    #out_path=Path("/home/appuser/venv/")
    root_folder = Path(sys.executable).parent.parent
    dctx=zstandard.ZstdDecompressor()
    with tempfile.TemporaryFile(suffix=".tar") as ofh:
        with filepath_final.open("rb") as ifh:
            dctx.copy_stream(ifh,ofh)
        ofh.seek(0)
        with tarfile.open(fileobj=ofh) as z:
            z.extractall(root_folder)
    os.system("ls /home/appuser/venv/lib/python3.9/lib-dynload")
    os.system("ldd /home/appuser/venv/lib/python3.9/lib-dynload/boost_python_meta_ext.so")
    dylib_folder = root_folder/f"lib/python{sys.version_info.major}.{sys.version_info.minor}/lib-dynload"
    sys.path.append("/home/appuser/venv/lib/python3.9/lib-dynload")
    sys.path.append("/home/appuser/venv/lib")
    os.system("rm /home/appuser/venv/lib/libstdc++.so.6")
    os.system("ln -s /home/appuser/venv/lib/libstdc++.so.6.0.30 /home/appuser/venv/lib/libstdc++.so.6")
    os.system("strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX")
    os.system("ldd /home/appuser/venv/lib/python3.9/lib-dynload/cctbx_xray_ext.so")
    os.system("strings /home/appuser/venv/lib/libstdc++.so.6 | grep GLIBCXX")
    #st.info(dylib_folder)
    os.environ["LD_LIBRARY_PATH"] = f"{dylib_folder.as_posix()}:{root_folder}/lib"
    #st.info(os.environ["LD_LIBRARY_PATH"])



#st.info(list(Path("/home/appuser/venv/").rglob("*tbx*")))
#import cctbx


#if Path("/home/appuser").exists():
#    # essential to avoid cctbx import errors
#    target = Path("/home/appuser/venv/share/cctbx")
#    if not target.exists():
#        target.symlink_to("/home/appuser/.conda/share/cctbx")
#
#    sys.path += ["/home/appuser/venv/lib/python3.9/lib-dynload"]
#    os.environ["PATH"] += os.pathsep + "/home/appuser/.conda/bin"


import numpy as np
import pandas as pd
from bokeh.plotting import ColumnDataSource, figure
from bokeh.models import Label, BasicTicker, ColorBar, LinearColorMapper, PrintfTickFormatter

from findmysequence_lib.findmysequence import fms_main
import mrcfile


def main():
    title = "map2seq: identification of proteins from density map"
    st.set_page_config(page_title=title, layout="wide")
    st.title(title)

    #https://discuss.streamlit.io/t/hide-titles-link/19783/4
    st.markdown(""" <style> .css-15zrgzn {display: none} </style> """, unsafe_allow_html=True)

    with st.expander(label="README", expanded=False):
            st.write("This is a Web App to help users identify proteins that best explain a density map. The user will provide a density map (mrc or ccp4 map format) and a main-chain model in pdb format using any sequence (e.g. Ala-model). The user can also specify a database of protein sequences such as human proteins or all proteins for the search. Currently, only the [findMySequence](https://journals.iucr.org/m/issues/2022/01/00/pw5018/) method is included although we plan to include additional methods in the future.  \nNOTE: the uploaded map/model files are **strictly confidential**. The developers of this app does not have access to the files")

    col1, [col2,col3] = st.sidebar, st.columns([1,0.5])

    mrc = None
    pdb = None

    with col1:

        mrc = None
        input_modes_map = {0:"upload", 1:"url", 2:"emd-xxxxx"}
        help_map = "Only maps in MRC (*\*.mrc*) or CCP4 (*\*.map*) format are supported. Compressed maps (*\*.gz*) will be automatically decompressed"
        input_mode_map = st.radio(label="How to obtain the input map:", options=list(input_modes_map.keys()), format_func=lambda i:input_modes_map[i], index=2, horizontal=True, help=help_map, key="input_mode_map")
        if input_mode_map == 0: # "upload a MRC file":
            label = "Upload a map in MRC or CCP4 format"
            help = None
            fileobj = st.file_uploader(label, type=['mrc', 'map', 'map.gz'], help=help, key="upload_map")
            if fileobj is not None:
                emd_id = extract_emd_id(fileobj.name)
                is_emd = emd_id is not None
                with open(os.path.join(tmpdir, fileobj.name), "wb") as f:
                    f.write(fileobj.getbuffer())
                mrc = tmpdir + "/" + fileobj.name
            else:
                return
        elif input_mode_map == 1: # "url":
            emd_id_default = "emd-10499"
            url_default = get_emdb_map_url(emd_id_default)
            help = "An online url (http:// or ftp://) or a local file path (/path/to/your/structure.mrc)"
            url = st.text_input(label="Input the url of a 3D map:", value=url_default, help=help, key="url_map").strip()
            if not url: return
            emd_id = extract_emd_id(url)
            is_emd = emd_id is not None and emd_id
            with st.spinner(f'Downloading {url}'):
                mrc = get_file_from_url(url)
            if mrc is None:
                st.warning(f"Failed to download [{url}]({url})")
                return
        elif input_mode_map == 2: # "emdb": randomly selects form emdb_ids_all
            with st.spinner(f'Downloading the list of all EMDB entries'):
                emdb_ids_all, resolutions = get_emdb_ids()
            if not emdb_ids_all:
                st.warning("failed to obtained a list of structures in EMDB")
                return
            url = "https://www.ebi.ac.uk/emdb/search/*%20?rows=10&sort=release_date%20desc"
            st.markdown(f'[All {len(emdb_ids_all):,} structures in EMDB]({url})')
            do_random_embid = st.checkbox("Choose a random EMDB ID", value=False, key="random_embid")
            if do_random_embid:
                help = "Randomly select another structure in EMDB"
                button_clicked = st.button(label="Change EMDB ID", help=help)
                if button_clicked:
                    import random
                    st.session_state.emd_id = 'emd-' + random.choice(emdb_ids_all)
            else:
                help = None
                label = "Input an EMDB ID (emd-xxxxx):"
                emd_id_default = "emd-23871"
                emd_id = st.text_input(label=label, value=emd_id_default, key='emd_id', help=help)
                if not emd_id: return
                emd_id = emd_id.lower().split("emd-")[-1]
                if emd_id not in emdb_ids_all:
                    import random
                    msg = f"EMD-{emd_id} is not a valid EMDB entry. Please input a valid id (for example, a randomly selected valid id 'emd-{random.choice(emdb_ids_all)}')"
                    st.warning(msg)
                    return
            if 'emd_id' in st.session_state: emd_id = st.session_state.emd_id
            else: emd_id = emd_id_default
            emd_id = emd_id.lower().split("emd-")[-1]
            url = get_emdb_map_url(emd_id)
            with st.spinner(f'Downloading EMD-{emd_id} from {url}'):
                mrc = get_file_from_url(url)
            if mrc is None:
                st.warning(f"Failed to download [EMD-{emd_id}]({url})")
                return
            resolution = resolutions[emdb_ids_all.index(emd_id)]
            msg = f'[EMD-{emd_id}](https://www.ebi.ac.uk/emdb/entry/EMD-{emd_id}) | resolution={resolution}Å'
            st.markdown(msg)

        if mrc is None or not Path(mrc).exists():
            st.warning(f"Failed to load density map")
            return

        mrc_changed = st.session_state.get("mrc_last", mrc) != mrc
        st.session_state.mrc_last = mrc

        st.markdown("""---""")

        #pdb input
        input_modes_model = {0:"upload", 1:"url", 2:"PDB ID"}
        help_model = "The input PDB model should have all backbone atoms (Cα,N,C) of each residue. Sidechain atoms are not required, resiudes can be labeled as any amino acids."
        input_mode_model = st.radio(label="How to obtain the input PDB file:", options=list(input_modes_model.keys()), format_func=lambda i:input_modes_model[i], index=2, horizontal=True, help=help_model, key="input_mode_model")
        pdb = None

        if input_mode_model == 0: # "upload a PDB file":
            label = "Upload a PDB file"
            fileobj = st.file_uploader(label, type=['pdb'], help=None, key="upload_model")
            if fileobj is not None:
                with open(os.path.join(tmpdir, fileobj.name), "wb") as f:
                    f.write(fileobj.getbuffer())
                pdb = tmpdir + "/" + fileobj.name
            else:
                return
        
        elif input_mode_model == 1: # "url":
            help = "An online url (http:// or ftp://) or a local file path (/path/to/your/model.pdb)"
            pdb_id_default = "6TGN"
            url_default = get_pdb_url(pdb_id_default)
            url = st.text_input(label="Input the url of a PDB model:", value=url_default, help=help, key="url_model").strip()
            if url:
                with st.spinner(f'Downloading {url}'):
                    pdb = get_file_from_url(url)
                    if pdb is None:
                        st.warning(f"Failed to download [{url}]({url})")
                        return
            else:
                return
        
        elif input_mode_model == 2: # "PDB ID":
            help = None
            label = "Input an PDB ID (for example: 4hhb):"
            pdb_id_default = "7mkf"
            pdb_id = st.text_input(label=label, key='pdb_id', value=pdb_id_default, help=help)
            pdb_id = pdb_id.upper()
            if pdb_id:
                pdb_url=get_pdb_url(pdb_id)
                with st.spinner(f'Downloading {pdb_id}.pdb from {pdb_url}'):
                    pdb = get_file_from_url(pdb_url)
                if pdb is None:
                    st.warning(f"Failed to download [PDB: {pdb_id}]({pdb_url})")
                    return
                msg = f'[PDB-{pdb_id}](https://www.rcsb.org/structure/{pdb_id})'
                st.markdown(msg)
            else:
                return

        if pdb is None or not Path(pdb).exists():
            st.warning(f"Failed to load the PDB model")
            return

        pdb_changed = st.session_state.get("pdb_last", pdb) != pdb
        st.session_state.pdb_last = pdb

        st.markdown("""---""")

        #db input
        input_modes_db = {0:"upload", 1:"url", 2:"human proteins", 3:"all curated proteins"}
        if not is_hosted(): input_modes_db[4] = "all proteins"
        help_db = "The input sequence database (.fa, .fa.gz, .fasta, or .fasta.gz)"
        input_mode_db = st.radio(label="Which sequence database to use:", options=list(input_modes_db.keys()), format_func=lambda i:input_modes_db[i], index=2, horizontal=True, help=help_db, key="input_mode_db")
        
        db = None
        info = "Searching {n:,} protein sequences"

        if input_mode_db == 0: # "upload":
            label = "Upload a fasta file (.fa, .fa.gz, .fasta, .fasta.gz)"
            if fileobj is None: return
            fileobj = st.file_uploader(label, type=['fa', 'fasta', 'fa.gz', 'fasta.gz'], help=None, key="upload_db")

            with open(os.path.join(tmpdir, fileobj.name), "wb") as f:
                f.write(fileobj.getbuffer())
            db = tmpdir + "/" + fileobj.name
        else:
            if input_mode_db == 1: # "url":
                help = "An online url (http:// or ftp://) or a local file path (/path/to/your/database.fa.gz)"
                url = st.text_input(label="Input the url of a sequence database (.fa, .fa.gz, .fasta, .fasta.gz):", help=help, key="url_db")
                if len(url)<1: return
            elif input_mode_db == 2: # "human proteins":
                url = "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
                info = "[{n:,} human protein sequences](https://www.uniprot.org/uniprotkb?facets=reviewed%3Atrue&query=%28proteome%3AUP000005640%29)"
            elif input_mode_db == 3: # "all curated proteins"
                url = "https://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz"
                info = "[{n:,} reviewed protein sequences](https://www.uniprot.org/uniprotkb?query=reviewed:true)"
            elif input_mode_db == 4: # "all proteins"
                #url = "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
                #info = "[{n:,} protein sequences (e.g. all known proteins) in Unreviewed (TrEMBL)](https://www.uniprot.org/help/downloads)"
                url = "https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/uniparc/uniparc_active.fasta.gz"
                info = "[{n:,} protein sequences (e.g. all known proteins) in UniParc](https://www.uniprot.org/help/uniparc)"
            with st.spinner(f'Downloading {url.strip()}'):
                db = get_file_from_url(url.strip())
        
        if db is None or not Path(db).exists():
            st.warning(f"Failed to load the protein sequence database")
            return      
        
        with st.spinner(f"Getting # of sequences in {db}"):
            n = number_of_sequences(db)
            st.markdown(info.format(n=n))

        direction_options = {0:"original", 1:"reversed"}
        help_direction=None
        direction_option = st.radio(label="Protein sequence direction:", options=list(direction_options.keys()), format_func=lambda i:direction_options[i], index=0, horizontal=True, help=help_direction, key="direction_option")
        
        handedness_options = {0:"original", 1:"flipped"}
        help_handedness=None
        handedness_option = st.radio(label="Map handedness:", options=list(handedness_options.keys()), format_func=lambda i:handedness_options[i], index=0, horizontal=True, help=help_handedness, key="handedness_option")
        
        if handedness_option in [1]: # flipped
            mrc, pdb = flip_map_model(mrc, pdb)

        if is_hosted():
            cpu = 1
        else:
            cpu = st.number_input("How many CPUs to use:", min_value=1, max_value=os.cpu_count(), value=2, step=1, help=f"a number in [1, {os.cpu_count()}]", key="cpu")

        st.markdown("""---""")
        run_button_clicked = st.button(label="Run")

        st.markdown("*Developed by the [Jiang Lab@Purdue University](https://jiang.bio.purdue.edu/map2seq). Report problems to Xiaoqi Zhang (zhang4377 at purdue.edu)*")

    if (mrc_changed or pdb_changed or input_mode_db in [3, 4]) and not run_button_clicked: return

    with col2:
        with st.spinner(info.format(n=number_of_sequences(db))):
            #remove_old_graph_log()
            seqin = None
            modelout = None
            res = map2seq_run(mrc, pdb, seqin, modelout, direction_option, handedness_option, db, cpu=cpu, outdir = tmpdir)
            if res is None:
                st.error(f"Failed")
                return

            xs, ys = res
                
        source = ColumnDataSource(data=dict(x=range(1,len(xs)+1),y=ys,ID=xs))
        top_source = ColumnDataSource(data=dict(x=[1],y=[ys[0]],ID=[xs[0]]))
        label = Label(x=1, y=ys[0], text=f'Best match: {xs[0]}', x_offset=10, y_offset=-5, text_font_size='16px', render_mode='canvas')

        TOOLTIPS = [('Rank','@x'),('Protein','@ID'),('E-val','@y')]
        p = figure(tooltips=TOOLTIPS, y_axis_type='log', title='')
        p.circle('x','y',source=source)
        p.circle('x','y',source=top_source, size=10,line_color='red',fill_color='red')
        p.yaxis.axis_label = 'E-values'
        p.xaxis.axis_label = 'Rank Order'
        p.y_range.flipped = True
        p.add_layout(label)
        p.xaxis.axis_label_text_font_size = "20pt"
        p.yaxis.axis_label_text_font_size = "20pt"
        p.xaxis.major_label_text_font_size = "16pt"
        p.yaxis.major_label_text_font_size = "16pt"
        
        st.bokeh_chart(p, use_container_width=True)
                
    with col2:
        df = pd.DataFrame({"E-val (log10)":np.log10(ys).T, "Protein":np.array(xs).T})
        df.index += 1
              
        def link_to_uniprot(s):
            pid = s.split('|')[1]
            url = f"https://www.uniprot.org/uniprotkb/{pid}"
            return f'<a href="{url}">{pid}</a>'            
        
        n = 10
        score_threshold = -2    # https://hmmer-web-docs.readthedocs.io/en/latest/searches.html
        n_df_col=4
        
        try:
            df.loc[:, "Uniprot ID"] = df.loc[:, "Protein"].str.split("|", expand=True).iloc[:, 1]
            df.loc[:, "URL"] = df.loc[:, "Protein"].apply(link_to_uniprot)
            df.loc[:, "Protein"] = df.loc[:, "Protein"].str.split("|", expand=True).iloc[:, -1]
            df.reset_index(inplace=True)
            df = df.rename(columns = {'index':'Rank'})
            df_top = df.iloc[:n, [0, 1, 2, 4]].copy()
        except:
            df.reset_index(inplace=True)
            df = df.rename(columns = {'index':'Rank'})
            df_top = df.iloc[:n, [0, 1, 2]].copy()
            n_df_col=3
        
        has_good_scores = df_top.iloc[:, 1].astype(float).min()<score_threshold
        def highlight_bad_score_rows(x, score_threshold=score_threshold):
            if x[1] > score_threshold:
                return ['background-color: red']*4
            else:
                return ['background-color: white']*4 
        df_top_style = df_top.style.apply(highlight_bad_score_rows, axis=1)
        df_top_style.hide(axis="index")
        st.markdown(f"**Top {n} matches:**")
        st.write(df_top_style.to_html(escape=False, index=False, justify="left"), unsafe_allow_html=True)
        if not has_good_scores:
            st.markdown(f":red[*No protein with meaningful score (<{score_threshold:g}). Check if the model is positioned properly in the map, or change the protein sequence database to **all curated proteins***]")

        st.download_button(
            label=f"Download the scores for {len(df):,} proteins",
            data=df.to_csv(index=False).encode('utf-8'),
            file_name='map2seq_results.csv',
            mime='text/csv'
        )

    with col2:
        tophit = xs[0]
        if st.button(label=f"Align top hit {tophit}", key="align_top_hit"):
            import pyfastx
            fa = pyfastx.Fasta(db)
            seqin = tmpdir+"/tmp.fasta"
            modelout = tmpdir+"/model_out.pdb"
            with open(tmpdir+"/tmp.fasta","w") as tmp:
                tmp.write(">"+xs[0]+"\n")
                tmp.write(fa[xs[0]].seq)
            
            with st.spinner("Processing..."):
                map2seq_run(mrc, pdb, seqin, modelout, direction_option, handedness_option, db, cpu=cpu, outdir = tmpdir)
            
            lines = []
            with open(tmpdir+"/seq_align_output.txt","r") as tmp:
                for line in tmp.readlines()[:-2]:
                    if "==>" in line or "Empty" in line or "WARNING" in line: continue
                    lines.append(line.rstrip())
            st.text("\n".join(lines))
                        
            with open(modelout,"r") as tmp:
                out_texts="".join(tmp.readlines())
                st.download_button("Download output model", data=out_texts, file_name=f"map2seq_model.pdb")
     
    with col2:
        with open(os.path.join(f'{tmpdir}/score_dict.pkl'),'rb') as f:
            score_dict_raw = pickle.load(f)
        
        st.download_button(
            label=f"Download the score matrix",
            data=score_dict_raw.to_csv().encode('utf-8'),
            file_name='score_matrix.csv',
            mime='text/csv'
        )        
            
        score_dict_raw.index.name="Residue"
        score_dict_raw.columns.name="AA"
        res_list=list(score_dict_raw.index)
        res_list.reverse()
        aa_list=list(score_dict_raw.columns)
            
        score_dict=pd.DataFrame(score_dict_raw.stack(),columns=["score"]).reset_index()
        
        
        
        colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
        mapper = LinearColorMapper(palette=colors, low=0, high=1)
        
        TOOLS = "hover,save,ywheel_pan,box_zoom,reset,wheel_zoom"
        
        hm = figure(title="Predicted Scores",
           #x_range=res_list, y_range=aa_list,
           #x_axis_location="below", width=900, height=400,
           x_range=aa_list, y_range=res_list,
           x_axis_location="above", width=900, height=len(res_list)*10,
           tools=TOOLS, toolbar_location='above',
           tooltips=[('Residue Position', '@Residue'), ('AA', '@AA'), ('Score','@score')])
        
        
        hm.grid.grid_line_color = None
        hm.axis.axis_line_color = None
        hm.axis.major_tick_line_color = None
        hm.axis.major_label_text_font_size = "7px"
        hm.axis.major_label_standoff = 0
        hm.xaxis.major_label_orientation = np.pi / 3
        #hm.rect(x="Residue", y="AA", width=1, height=1,
        #   source=score_dict,
        #   fill_color={'field': 'score','transform': mapper},
        #   line_color=None)
        
        hm.rect(x="AA", y="Residue", width=1, height=1,
           source=score_dict,
           fill_color={'field': 'score','transform': mapper},
           line_color=None)
        
        color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="7px",
                     ticker=BasicTicker(desired_num_ticks=len(colors)),
                     formatter=PrintfTickFormatter(format="%.1f"),
                     label_standoff=6, border_line_color=None)
        hm.add_layout(color_bar, 'right')
        
        hm
        
    
    with col3:
        #st.markdown("PDB Model Overview")
        #plot_pdb_model(pdb)
        
        st.markdown("Map Projections")
        plot_density_projection(mrc)

    remove_old_pdbs(keep=10)
    remove_old_maps(keep=10)
    #remove_old_graph_log()

#@st.experimental_memo(persist='disk', max_entries=1, ttl=60*60*24, show_spinner=False, suppress_st_warning=True)
#def plot_pdb_model(pdb):
#    from stmol import showmol
#    import py3Dmol
#    with open(pdb) as f:
#        s="".join([line for line in f])
#    s_view=py3Dmol.view(data=s,width=400,height=300)
#    s_view.setStyle({'cartoon':{'color':'spectrum'}})
#    showmol(s_view)    

#@st.experimental_memo(persist='disk', max_entries=1, ttl=60*60*24, show_spinner=False, suppress_st_warning=True)
def plot_density_projection(mrc):
    mrc_data = mrcfile.open(mrc, 'r+')
    v_size=mrc_data.voxel_size
    nx=mrc_data.header['nx']
    ny=mrc_data.header['ny']
    nz=mrc_data.header['nz']
    apix=v_size['z']
    
    data=mrc_data.data
    ret = data.sum(axis=0)
    ret=normalize(ret)
    st.image(ret)

    ret = data.sum(axis=1)
    ret=normalize(ret)
    st.image(ret)    

    ret = data.sum(axis=2)
    ret=normalize(ret)
    st.image(ret)
    
    #mrc_data = mrcfile.open(mrc, 'r+')
    #v_size=mrc_data.voxel_size
    #nx=mrc_data.header['nx']
    #ny=mrc_data.header['ny']
    #nz=mrc_data.header['nz']
    #apix=v_size['z']
    
    #data=mrc_data.data
        
    #X,Y,Z=np.mgrid[0:apix*nx:nx*1j,0:apix*ny:ny*1j,0:apix*nz:nz*1j]

        
    #import plotly.graph_objects as go
    ##mrc_fig=go.Figure(data=go.Volume(x=np.arange(0,nx*apix,apix),y=np.arange(0,ny*apix,apix),z=np.arange(0,nz*apix,apix),value=data,isomin=0.1,isomax=0.8,opacity=1,surface_count=200))
    #mrc_fig=go.Figure(data=go.Volume(x=X.flatten(),y=Y.flatten(),z=Z.flatten(),value=data.flatten(),isomin=0.1,isomax=0.8,opacity=0.1,surface_count=20))
    #st.plotly_chart(mrc_fig,use_container_width=True)

def normalize(data, percentile=(0, 100)):
    p0, p1 = percentile
    vmin, vmax = sorted(np.percentile(data, (p0, p1)))
    data2 = (data-vmin)/(vmax-vmin)
    return data2
       
def remove_old_graph_log():
    dir = os.listdir(tmpdir)
    for item in dir:
        if item.endswith(".png") or item.endswith(".txt"):
            os.remove(os.path.join(tmpdir, item))

def remove_old_maps(keep=0):
    map_files = [os.path.join(tmpdir, item) for item in os.listdir(tmpdir) if item.endswith(".mrc") or item.endswith(".map") or item.endswith(".map.gz")]
    if keep>0:
        map_files = sorted(map_files, key=lambda f: os.path.getmtime(f))[:-keep]
    for f in map_files:
        os.remove(os.path.join(tmpdir, f))

@st.experimental_memo(persist=True, show_spinner=False)
def number_of_sequences(db_fasta):
    import pyfastx
    fa = pyfastx.Fasta(db_fasta)
    return len(fa)

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

# do not use st cache
def get_file_from_url(url):
    local_file_name = Path(url).name
    local_file_name = Path(tmpdir)/local_file_name
    if local_file_name.suffix == ".gz":
        filename_final = local_file_name.parent / local_file_name.stem
    else:
        filename_final = local_file_name    
    if filename_final.exists():
        return filename_final.as_posix()
    
    if is_jianglab():
        db_folder = Path(sys.executable).parent.parent.parent.parent / "protein_sequence_db"
        db_file = db_folder / filename_final.name
        if db_file.exists():
            filename_final.symlink_to(db_file)
            return filename_final.as_posix()
        db_file = db_folder / local_file_name.name
        if db_file.exists():
            local_file_name.symlink_to(db_file)
    
    if not local_file_name.exists():
        url_final = get_direct_url(url)    # convert cloud drive indirect url to direct url
        ds = np.DataSource(None)
        if not ds.exists(url_final):
            st.error(f"ERROR: {url} could not be downloaded. If this url points to a cloud drive file, make sure the link is a direct download link instead of a link for preview")
            st.stop()
        with ds.open(url_final) as fp:
            local_file_name = Path(tmpdir)/Path(fp.name).name
            import shutil
            shutil.move(fp.name, local_file_name)

    if local_file_name.suffix == ".gz":
        import gzip, shutil
        with gzip.open(local_file_name, 'r') as f_in, open(filename_final, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return filename_final.as_posix()

def extract_emd_id(text):
    import re
    pattern = '.*emd_([0-9]*)\.map.*'

    match = re.search(pattern, text, re.IGNORECASE)
    if match:
        emd_id = match.group(1)
    else:
        emd_id = None
    return emd_id

@st.experimental_memo(max_entries=1, ttl=60*60*24, show_spinner=False, suppress_st_warning=True)
def get_emdb_ids():
    try:
        import pandas as pd
        entries = pd.read_csv("https://www.ebi.ac.uk/emdb/api/search/current_status:%22REL%22%20?wt=csv&download=true&fl=emdb_id,resolution")
        emdb_ids = list(entries.iloc[:,0].str.split('-', expand=True).iloc[:, 1].values)
        resolutions = entries.iloc[:,1].values
    except:
        emdb_ids = []
        resolutions = []
    return emdb_ids, resolutions

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
	
#-------------------------------End Map Functions-------------------------------

#-------------------------------Model Functions-------------------------------
def remove_old_pdbs(keep=0):
    pdb_files = [os.path.join(tmpdir, item) for item in os.listdir(tmpdir) if item.endswith(".pdb") or item.endswith(".pdb.gz")]
    if keep>0:
        pdb_files = sorted(pdb_files, key=lambda f: os.path.getmtime(f))[:-keep]
    for f in pdb_files:
        os.remove(os.path.join(tmpdir, f))

@st.experimental_memo(max_entries=1, ttl=60*60*24*7, show_spinner=False, suppress_st_warning=True)
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

def flip_map_model(map_name, pdb_name):
    map_flip = Path(map_name).with_suffix(".flip.mrc")
    if not map_flip.exists():
        with mrcfile.open(map_name) as mrc_data:
            apix=mrc_data.voxel_size.z
            with mrcfile.new(str(map_flip), overwrite=True) as mrc_data_flip:
                mrc_data_flip.set_data(mrc_data.data[::-1, :, :])   # z-flip
                mrc_data_flip.voxel_size = apix

    pdb_flip = Path(pdb_name).with_suffix(".flip.pdb")
    if pdb_flip.exists():
        return str(map_flip), str(pdb_flip)

    with mrcfile.open(map_name, header_only=True) as mrc_data:
        apix=mrc_data.voxel_size.z
        nz=mrc_data.header['nz']

    with open(pdb_name,'r') as f:
        lines=f.readlines()
        orig_atom_lines = []
        for line in lines:
            if line[0:4] == "ATOM":
                orig_atom_lines.append(line)
            elif line[0:3] == "TER":
               orig_atom_lines.append(line)
    with open(pdb_flip, 'w') as o:
        for line in orig_atom_lines:
            if line[0:3] != "TER":
                coord_vec = np.array(line[30:54].split()).astype(np.float32)
                new_z=(nz-1)*apix-coord_vec[2]
                line = list(line)
                line[30:54] = list(" " + f'{coord_vec[0]:>7.3f}' + " " + f'{coord_vec[1]:>7.3f}' + " " + f'{new_z:>7.3f}')
                line = "".join(line)
            o.write(line)
    return str(map_flip), str(pdb_flip)

@st.experimental_memo(max_entries=10, ttl=60*60, show_spinner=False, suppress_st_warning=True)
def map2seq_run(map, pdb, seqin, modelout, rev, flip, db, cpu=2, outdir="tempDir/"):
    os.environ['cpu'] = f"{cpu}"

    map = os.path.abspath(map)
    pdb = os.path.abspath(pdb)
    db = os.path.abspath(db)
    outdir = os.path.abspath(outdir)

    if outdir[-1] != "/":
        outdir += "/"

    basename = f"map2seq_fms"

    fms_main.fms_run(mapin=map, modelin=pdb, seqin=seqin, modelout=modelout, db=db, tmpdir=outdir, outdir=outdir, rev=rev, flip=flip, tophits=np.iinfo(np.uint32).max)
    
    hmmer_out = "hmmer_output.txt"
    num = parse_file(f"{outdir}{basename}.png", f"{outdir}{hmmer_out}")
    if num == -1:
        return None # failed

    with open(os.path.join(outdir, f'{basename}.png_x.pkl'),'rb') as inf:
        xs = pickle.load(inf)
    with open(os.path.join(outdir, f'{basename}.png_y.pkl'),'rb') as inf:
        ys = pickle.load(inf)

    return (xs, ys)
        
def make_graph(ids, e_vals, outputFile):
    
    #output_file('{}.html'.format(outputFile))

    #source = ColumnDataSource(data=dict(x=range(len(ids)),y=e_vals,ID=ids))
    #top_source = ColumnDataSource(data=dict(x=[0],y=[e_vals[0]],ID=[ids[0]]))
    #label = Label(x=0, y=e_vals[0], text='Best Match', x_offset=10, y_offset=-5, render_mode='canvas')
  
    #TOOLTIPS = [('index','$index'),('ID','@ID'),('E-val','@y')]
   
    #p = figure(width=400,height=400,tooltips=TOOLTIPS,y_axis_type='log', title='Ranked Sequences')
    #p.circle('x','y',source=source)
    #p.circle('x','y',source=top_source, size=10,line_color='red',fill_color='red')
    #p.yaxis.axis_label = 'E-values'
    #p.xaxis.axis_label = 'Rank Order'
    #p.y_range.flipped = True
    #p.add_layout(label)
    
    #save(p)
    
    with open('{}_x.pkl'.format(outputFile),'wb') as o:
        pickle.dump(ids,o,pickle.HIGHEST_PROTOCOL)
    with open('{}_y.pkl'.format(outputFile),'wb') as o:
        pickle.dump(e_vals,o,pickle.HIGHEST_PROTOCOL)


def parse_file(outputFile, filepath):
    matches_found = "Best matches"
    no_matches_found = "No matches found"

    ids = []
    e_vals = []
    with open(filepath, 'r') as file:
        firstline = file.readline().rstrip()
        if firstline.find(matches_found) == -1:
            print(no_matches_found)
            return -1
        for line in file:
            line = line.rstrip()
            #list = line.split('|')
            ##list[0:3] = ["|".join(list[0:3])]
            ###list = line.split(' ')
            ##list[0] = list[0].strip()
            ##list[1] = list[1].removeprefix('E-value=')
            ##list[1] = float(list[1])
            #curr_id="|".join(list[:-1])
            #curr_ev=float(list[-1].removeprefix('E-value='))
            list=line.split()
            curr_id=list[0]
            curr_ev=float(list[-1])
            ids.append(curr_id)
            e_vals.append(curr_ev)
                        
        make_graph(ids, e_vals, outputFile)
    return 1

def is_hosted():
    ret = Path("/home/appuser").exists()
    return ret

def is_jianglab():
    ret = Path("/net/jiang").exists()
    return ret

@st.cache(persist=True, show_spinner=False)
def setup_anonymous_usage_report():
    try:
        import pathlib, stat
        index_file = pathlib.Path(st.__file__).parent / "static/index.html"
        index_file.chmod(stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IROTH)
        txt = index_file.read_text()
        if txt.find("gtag/js?")==-1:
            txt = txt.replace("<head>", '''<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-VSTDDFT4HW"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag('js', new Date());gtag('config', 'G-VSTDDFT4HW');</script>''')
            index_file.write_text(txt)
    except:
        pass


if __name__ == "__main__":
    setup_anonymous_usage_report()

    if is_hosted():
        ## essential to avoid cctbx import errors
        #target = Path("/home/appuser/venv/share/cctbx")
        #if not target.exists():
        #    target.symlink_to("/home/appuser/.conda/share/cctbx")

        sys.path += ["/home/appuser/venv/lib/python3.9/lib-dynload"]
        #os.environ["PATH"] += os.pathsep + "/home/appuser/.conda/bin"

    main()
