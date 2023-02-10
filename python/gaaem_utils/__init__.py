import numpy as np;
import pandas as pd
import sys
import matplotlib.pyplot as plt
import copy, os


def calc_lineOffset(data):
    if ('UTMX' in data['flightlines'].keys()) and ('UTMY' in data['flightlines'].keys()):
        pos_keys=['UTMX', 'UTMY']
    elif ('utmx' in data['flightlines'].keys()) and ('utmy' in data['flightlines'].keys()):
        pos_keys=['utmx', 'utmy']
    elif ('easting' in data['flightlines'].keys()) and ('northing' in data['flightlines'].keys()):
        pos_keys=['easting', 'northing']
    else:
        raise Exception("Sorry, no coordinates with labels UTMX, UTMY or utmx, utmy, found in data") 
    dx=data['flightlines'][pos_keys[0]].diff()
    dx.iloc[0]=0
    dy=data['flightlines'][pos_keys[1]].diff()
    dy.iloc[0]=0
    data['flightlines'].insert(len(data['flightlines'].columns), 'lineoffset', np.cumsum( np.sqrt(dx**2 + dy**2) ) )

def make_XYZ_model(xyz):
    for key in xyz["layer_data"].keys():
            clear_column_name(xyz["layer_data"][key])
    xyz["layer_data"]["rho_i"]=1./xyz["layer_data"]["conductivity"]
    xyz["layer_data"]["dep_bot"]=xyz["layer_data"]['thickness'].cumsum(axis=1)
    xyz["layer_data"]["dep_top"]=xyz["layer_data"]["dep_bot"]-xyz["layer_data"]['thickness']
    xyz["flightlines"]["alt"]=xyz["flightlines"].tx_height  
    if "inverted_tx_height" in xyz["flightlines"].columns:
        xyz["flightlines"]["invalt"]=xyz["flightlines"].inverted_tx_height
    else:
        xyz["flightlines"]["invalt"]=xyz["flightlines"].tx_height * np.nan


def plot_model_section(model, ax, keyx="utmy", res_key='rho_i', elev_key='elevation', doi_key='doi_standard', 
                       cmap="jet", clim=[-1, 4], hideBelowDOI=True, attr='log_res', cb_orientation='vertical', showRMS=False):
    # some preparations
    if ('lineoffset' in keyx) and not('lineoffset' in model['flightlines'].columns):
        calc_lineOffset(model)
    
    local_x = model['flightlines'][keyx]
    if attr=='log_res':
        image_data=np.log10(model['layer_data'][res_key].values.T)
        cb_label=r'Resistivity [log10($\Omega$m)]'
    elif attr=='res':
        image_data=model['layer_data'][res_key].values.T
        cb_label=r'Resistivity [$\Omega$m]'
    elif attr=='cond':
        image_data=1./model['layer_data'][res_key].values.T
        cb_label=r'Conductivity [S/m]'
    elif attr=='log_cond':
        image_data=np.log10(1./model['layer_data'][res_key].values.T)
        cb_label=r'Conductivity [log10(S/m)]'
    elif attr=='log_cond_mSm':
        image_data=np.log10(1./model['layer_data'][res_key].values.T*1000)
        cb_label=r'Conductivity [log10(mS/m)]'
    depth_top=model['layer_data']['dep_top'].values.T
    if model['layer_data']['dep_top'].shape[1]==model['layer_data']['dep_bot'].shape[1]+1:
        df=model['layer_data']['dep_bot']
        ncol=df.shape[1]
        maxDepth=df.loc[:,ncol-1]+50
        model['layer_data']['dep_bot'].insert(ncol, df.columns[-1]+1, maxDepth)

    depth_bot=model['layer_data']['dep_bot'].values.T
    depth=model['flightlines'][elev_key].values - (depth_top+depth_bot)/2
    # actual plotting
    pm=ax.pcolormesh(local_x, depth, image_data, cmap=cmap, shading='auto', vmin=clim[0], vmax=clim[1] )
    ax.set_ylabel('Elevation (m)')
    ax.set_xlabel(keyx)
    cb=plt.colorbar(pm, ax=ax, orientation=cb_orientation, shrink=0.5 )
    cb.set_label(cb_label)
    ax.grid(linestyle='--', linewidth=0.5)
    ax.plot(local_x, model['flightlines'][elev_key] + model['flightlines'].alt, 'b-', label="alt.", lw=.5)
    ax.plot(local_x, model['flightlines'][elev_key] + model['flightlines'].invalt, 'g--', label="inv. alt.", lw=.5)
    #ax.plot(local_x, model['flightlines'][elev_key] - model['flightlines'][doi_key], "k-", label="DOI", lw=.5)
    ax.legend()
    if hideBelowDOI:
        z=(model['flightlines'][elev_key] - model['flightlines'][doi_key]).values
        #zmax=(model['flightlines'][doi_key] - model['layer_data']['dep_bot'].max(axis=1)).values-100
        zmax=(model['flightlines'][elev_key] - model['layer_data']['dep_bot'].max(axis=1)).values-50
        poly_coords=[]
        xbuffer=0
        for n, x in enumerate(local_x):
            if n==1:
                poly_coords.append((x+xbuffer, z[n]))
            elif n==len(local_x)-1:
                poly_coords.append((x-xbuffer, z[n]))
            else:
                poly_coords.append((x, z[n]))
        
        if local_x.iloc[0]<local_x.iloc[-1]:
            poly_coords.append( (local_x.max(), np.max(zmax) ) )
            poly_coords.append( (local_x.min(), np.max(zmax) ) )
        else:
            poly_coords.append( (local_x.min(), np.max(zmax) ) )
            poly_coords.append( (local_x.max(), np.max(zmax) ) )
        ax.add_patch(plt.Polygon(poly_coords, color='white', alpha=0.5) )
    if showRMS:
        ax2=ax.twinx()
        ax2.plot(local_x, model['flightlines'].resdata)
        ax2.set_ylabel('data RMS')
        ax2.set_ylim([-6, 3])
    ax.set_xlim([local_x.min(), local_x.max()])


def clear_column_name(df):
    rename_dict={}
#    for c in df.columns:
#        rename_dict[c]=c.split("_")[-1]
    for n, c in enumerate(df.columns):
        rename_dict[c]=str(n)    
    df.rename(columns=rename_dict, inplace=True)

    
def read_GAAEM_data(data_fullfile):
    if ".dat" in data_fullfile:
        header_fullfile=data_fullfile.split(".dat")[0]+".hdr"
    elif ".asc" in data_fullfile:
        header_fullfile=data_fullfile.split(".asc")[0]+".hdr"
    with open(header_fullfile) as fid:
        lines=fid.readlines()
    colnames={}
    multicolnames={}
    for line in lines:
        col, name= line.split()
        if len(col.split("-"))==1:
            colnames[int(col)-1]=name
        elif len(col.split("-"))==2:
            col_beg, col_end = col.split("-")
            multicolnames[name]=[]
            for c in np.arange(int(col_beg), int(col_end)+1):
                colname="_".join([name, str(c-int(col_beg))])
                colnames[c-1]=colname
                multicolnames[name].append(colname)
                
                    
                
    df=pd.read_csv(data_fullfile, header=None, sep="\s+").rename(columns=colnames)
    layer_data={}
    for key in multicolnames.keys():
        layer_data[key]=pd.concat([df.pop(x) for x in multicolnames[key]], axis=1)
        #clear_column_name(layer_data[key])
        
    data={"flightlines":df,
          "layer_data":layer_data
         }
        
    return data


def get_GALEI_files(data_dir, data_prefix="inversion.output", data_suffix=".asc"):
    file_list=[]
    for root, dirs, files in os.walk(data_dir):
            for file in files:
                if file.endswith(data_suffix) and file.startswith(data_prefix):
                    file_list.append(os.path.join(root, file))
    if len(file_list)==0:
        raise Exception("No inversion files found in this folder")
    return file_list

def read_GALEI_files(data_dir, data_prefix="inversion.output", data_suffix=".asc"):
    file_list = get_GALEI_files(data_dir, data_prefix="inversion.output", data_suffix=".asc")
    print("Reading files:{}".format(file_list))
    alldata={}
    for file in file_list:
        data=read_GAAEM_data(file)
        if len(alldata.keys())==0:
            alldata=data
        else:
            alldata["flightlines"]=pd.concat([alldata["flightlines"], data["flightlines"]], axis=0)
            for key in data["layer_data"].keys():
                if key in alldata["layer_data"].keys():
                    alldata["layer_data"][key] = pd.concat([alldata["layer_data"][key], data["layer_data"][key]], axis=0)
                else:
                    print("something is wrong")
                    print("cannot find key: {0} in file: {1}".format(key,file))
    df=copy.deepcopy(alldata["flightlines"]) 
    layer_data_columns={}
    for key in alldata["layer_data"].keys():
        df=pd.concat([df, alldata["layer_data"][key]], axis=1)
        layer_data_columns[key]=alldata["layer_data"][key].columns

    df.sort_values("uniqueid", inplace=True)
    df.reset_index(inplace=True, drop=True)

    layer_data={}
    for key in layer_data_columns.keys():
        layer_data[key]=pd.concat([df.pop(x) for x in layer_data_columns[key]], axis=1)
        #clear_column_name(layer_data[key])

    data_out={"flightlines":df,
          "layer_data":layer_data
         }
    return data_out
    
def GALEI_invQCplot(model, figsize=(8,11), keyx="lineoffset"):
    if (keyx =="lineoffset") and not("lineoffset" in model["flightlines"].columns):
        calc_lineOffset(model)
    elif  not(keyx in model["flightlines"].columns):
        raise Exception("keyx not found in model['flightlines']")
    fig, ax = plt.subplots(4,1,figsize=figsize, sharex=True)
    ax[0].plot(model["flightlines"][keyx], 
              model["layer_data"]["observed_EMSystem_1_ZS"].abs(),
              ".", ms=2,
              label="observed_EMSystem_1_ZS")
    ax[0].set_prop_cycle(None)
    ax[0].plot(model["flightlines"][keyx], 
              model["layer_data"]["predicted_EMSystem_1_ZS"].abs(),
              "-", lw=1,
              label="predicted_EMSystem_1_ZS")
    ax[0].set_yscale("log")
    ax[0].set_title("data LM")

    ax[1].plot(model["flightlines"][keyx], 
              model["layer_data"]["observed_EMSystem_2_ZS"].abs(),
               ".", ms=2,
              label="observed_EMSystem_2_ZS")
    ax[1].set_prop_cycle(None)
    ax[1].plot(model["flightlines"][keyx], 
              model["layer_data"]["predicted_EMSystem_2_ZS"].abs(),
               "-", lw=1,
              label="predicted_EMSystem_2_ZS")
    ax[1].set_yscale("log")
    ax[1].set_title("data HM")

    plot_model_section(model, ax=ax[2], keyx=keyx, hideBelowDOI=False, cb_orientation="horizontal", clim=[0, 3.5])
    ax[1].set_title("inverted model")

    for phi  in ["PhiM", "PhiD", "PhiC", "PhiT", "PhiG", "PhiS"]:
        ax[3].plot(model["flightlines"][keyx], model["flightlines"][phi], label=phi)
    ax[3].set_yscale("log")
    ax[3].set_xlabel(keyx)
    ax[3].grid()
    ax[3].legend()

    plt.tight_layout()
    return fig, ax