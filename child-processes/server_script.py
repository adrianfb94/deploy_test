import numpy as np
import sys
import subprocess
import datetime
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning)
filterwarnings(action='ignore', category=UserWarning)
# from netCDF4 import Dataset, num2date
import struct
import os, glob
import xarray as xr
import numpy.ma as ma

# pip install python-dotenv
from dotenv import load_dotenv

# pip install gitpython
# import git


load_dotenv()


# os.environ['GADDIR'] = '/usr/lib/cgi-bin/grads-2.1.a2.oga.1/Classic/data'
# os.environ['GASCRP'] = '/usr/lib/cgi-bin/grads-2.1.a2.oga.1/Classic/scripts'
# os.environ['GAUDFT'] = '/usr/lib/cgi-bin/grads-2.1.a2.oga.1/Classic/data'


my_id = os.getenv('ID')
username = os.getenv('USERNAME')
password = os.getenv('PASSWORD')
token = os.getenv('TOKEN')
home_server = os.getcwd()
repo_name = 'deploy_test'


##print('HOME is: {}'.format(home_server))
##print('repo_name is: {}'.format(repo_name))
##print("Mi numero de carne es: {}".format(my_id))
##print("Mi username: {}".format(username))
##print("Mi password: {}".format(password))
##print("Mi token: {}".format(token))


# os.system('git config user.email "adrianfuentesbarrios@gmail.com"')
# os.system('git config user.name "adrianfb94"')

# url = f"https://{token}@github.com/{username}/{repo_name}.git"

# #print('tengo el REPO')



# file = 'another_file_added_from_server.txt'
# if not os.path.exists('/'.join([home,file])):
#     ##print("no existe {}".format('/'.join([home,file])))

#     subprocess.run('touch {}'.format('/'.join([home,file])), shell=True)
#     ##print('{} was created'.format(file))

# else:
#     ##print("si existe {}".format('/'.join([home,file])))

#     ##print('estoy haciendo el init')
#     subprocess.run('git init', shell=True, capture_output=True)
#     # os.system('git init')
#     ##print('ya hice el init')

#     ##print("estoy haciendo el add")
#     subprocess.run('git add {}'.format(file), shell=True, capture_output=True)
#     # os.system('git add file_added_from_server.txt')
#     ##print('ya hice el add')

#     ##print("estoy haciendo el commit")
#     subprocess.run('git commit -m "added {} from SERVER"'.format(file), shell=True, capture_output=True)
#     # os.system('git commit -m "added from SERVER"')
#     ##print('ya hice el commit')

#     ##print("estoy haciendo el push")
#     subprocess.run('git push https://{}@github.com/{}/{}.git HEAD:main'.format(token, username, repo_name), shell=True, capture_output=True)
#     # os.system(f'git push https://{token}@github.com/{username}/{repo_name}.git HEAD:main')
#     ##print('ya hice el push')











fechahora = str(datetime.datetime.now())
fecha = fechahora.split(' ')[0].replace('-','')
hora = fechahora.split(' ')[1].split('.')[0].replace(':','_')
nowdate = '_'.join([fecha, hora])

# def addpng(cad):
#     if (len(cad.split('.'))==1):
#         cad += '.png'
#     return cad

# import numpy.ma as ma
# def get_nc(var,f, miss):
#     nc_file = Dataset(f)
#     nc_lat = nc_file.variables['lat'][:]
#     nc_lon = nc_file.variables['lon'][:]
#     nc_var = nc_file.variables[var][:]
#     nc_var = nc_file.variables[var][:]
#     nc_var.fill_value = miss
#     nc_var[nc_var.mask] = nc_var.fill_value


#     nc_time_array = np.asarray(nc_file.variables['time'])

#     nc_time_units = nc_file.variables['time'].units
#     nc_file.close()
#     return nc_lat, nc_lon, nc_var, nc_time_array, nc_time_units


def get_nc(var,f, miss):

    # file = subprocess.run('tar -xzvf {}'.format(f), shell=True, capture_output=True).stdout.decode('utf-8').split('\n')[0]

    # nc_file = Dataset(f)
    # nc_lat = nc_file.variables['lat'][:]
    # nc_lon = nc_file.variables['lon'][:]
    # if nc_lon.min()>0 and nc_lon.max()>0:
    #     nc_lon = nc_lon-360

    # nc_var = nc_file.variables[var][:]
    # nc_var.fill_value = miss
    # nc_var[nc_var.mask] = nc_var.fill_value


    # nc_time_array = np.asarray(nc_file.variables['time'])

    # nc_time_units = nc_file.variables['time'].units
    # nc_file.close()
    # return nc_lat, nc_lon, nc_var, nc_time_array, nc_time_units
    ds = xr.load_dataset(f)
    # ds[var] = ds.variables[var].fillna(miss)
    # ##print(ds[var])
    # exit()

    if ds.lon.values.min()>0 and ds.lon.values.max()>0:
        ds['lon'] = ds['lon']-360

    values = ds[var].fillna(miss)
    # values = ds[var].fillna(miss)
    times = ds.time.values

    subprocess.run('rm child-processes/data/data4drought/data/echam5/*.nc', shell=True, capture_output=False)

    return ds, ds.lat, ds.lon, values, times
    # return ds, ds.lat, ds.lon, ds[var], times



def get_ip(raw_ip):
    if len(raw_ip.split(':'))==1:
        return raw_ip
    else:
        return raw_ip.split(':')[-1]


# cad1 = 'latiN=29&longW=-94.6&longE=-50&latiS=5&modelo=3&indice=1&periodoS=3&annoIRef=1961&annoFRef=1980&annoIEva=1961&annoFEva=1980&moddro=-1&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=4&comportamiento=0&lsm=0'
# cad2 = 'latiN=29&longW=-94.6&longE=-50&latiS=5&modelo=3&indice=1&periodoS=2&annoIRef=1961&annoFRef=1980&annoIEva=1961&annoFEva=1980&moddro=-1.0&sevdro=-1.5&extdrou=-2.0&umbral=2&def=1&comportamiento=0&lsm=0'
# cad3 = 'latiN=29&longW=-94.5&longE=-50&latiS=5&modelo=3&indice=1&periodoS=3&annoIRef=1961&annoFRef=2000&annoIEva=1976&annoFEva=1990&moddro=-1&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=4&comportamiento=0&lsm=1'
# cad4 = 'latiN=8&longW=-70&longE=-70&latiS=8&modelo=3&indice=1&periodoS=1&annoIRef=1961&annoFRef=1978&annoIEva=1961&annoFEva=2002&moddro=-1&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=1&comportamiento=0&lsm=1'

# # # Arnoldo first line
# line1="latiN=29&longW=-94.6&longE=-50&latiS=5&modelo=3&indice=1&periodoS=1&annoIRef=1961&annoFRef=1980&annoIEva=1961&annoFEva=1980&moddro=-1&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=2&comportamiento=0&lsm=0"

# # Arnoldo second line
# line2="latiN=29&longW=-94.5&longE=-50&latiS=5&modelo=3&indice=1&periodoS=3&annoIRef=1961&annoFRef=2000&annoIEva=1976&annoFEva=1990&moddro=-1&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=4&comportamiento=1&lsm=0"


# # # # Arnoldo third line
# line3="latiN=29&longW=-94.6&longE=-50&latiS=5&modelo=3&indice=2&periodoS=1&annoIRef=1961&annoFRef=1980&annoIEva=1961&annoFEva=1980&moddro=-1&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=2&comportamiento=0&lsm=0"



# *******Aqui empieza la fiesta*********

# cad = 'latiN=29&longW=-94.6&longE=-50&latiS=5&modelo=3&indice=2&periodoS=1&annoIRef=1961&annoFRef=1980&annoIEva=1961&annoFEva=1980&moddro=-1&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=2&comportamiento=0&lsm=0'
# argumentos = cad.split('&')
# ip = '127.0.0.1'

# cgiline = 'latiN=29&longW=-94.6&longE=-50&latiS=5&modelo=3&indice=1&periodoS=1&annoIRef=1961&annoFRef=1980&annoIEva=1961&annoFEva=1980&moddro=-1&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=2&comportamiento=0&lsm=0'
# argumentos = cgiline.split('&')
# ip = '127.0.0.1'






# SERVER
server_cad = sys.argv[1:][0]
argumentos = server_cad.split('&')
client_ip = sys.argv[1:][1]
ip = get_ip(client_ip)


# # LOCAL
# server_cad = 'latiN=29&longW=-94.6&longE=-50&latiS=5&modelo=3&indice=1&periodoS=1&annoIRef=1961&annoFRef=1980&annoIEva=1961&annoFEva=1980&moddro=-1&sevdro=-1.5&extdrou=-2&umbral=-0.5&def=2&comportamiento=0&lsm=1'
# argumentos = server_cad.split('&')
# ip = '127.0.0.1'

dict = {}
for item in argumentos:
    [name, value] = item.split('=')
    dict[name] = value

# dict = {}
# for i, arg in enumerate(argumentos):
#     dict[listvar[i]]=arg.split('=')[1]


# for i in range(len(listvar)):
#     ##print(listvar[i], dict[listvar[i]])
# exit()

models = ['INIT_0',
          'crudata', 
          'udeldata', 
          'echam5', 
          'aenwh', 
          'aexsa', 
          'aexsc', 
          'aexsk', 
          'aexsl', 
          'aexsm', 
          'ensmod', 
          'rcp26', 
          'rcp45', 
          'rcp85'
          ]



# listvar = ['latiN', 
#            'longW',
#            'longE',
#            'latiS',
#            'modelo',
#            'indice',
#            'periodoS',
#            'annoIRef',
#            'annoFRef',
#            'annoIEva',
#            'annoFEva',
#            'moddro',
#            'sevdro',
#            'extdrou',
#            'umbral',
#            'def',
#            'comportamiento',
#            'lsm'
#            ]

dindex = int(dict['indice']) #nombre del file de salida para ser procesado
modeldir = int(dict['modelo']) #nombre del file de entrada para ser preparado
modeldir2 = models[modeldir]
ilon = float(dict['longW']) #longitud inicial
elon = float(dict['longE']) #longitud final (no puede ser la misma que ilon)
ilat = float(dict['latiS']) #latitud inicial
elat = float(dict['latiN']) #latitud final (no puede ser la misma que ilat)
iryear = int(dict['annoIRef']) #anno inicial referencia
eryear = int(dict['annoFRef']) #anno final referencia
iayear = int(dict['annoIEva']) #anno inicial assessment
eayear = int(dict['annoFEva']) #anno final assessment
avg = int(dict['comportamiento']) # 0 false (no avg) , 1 true (si avg)
spilen = int(dict['periodoS']) #spi lenght
spilim = float(dict['umbral']) # spi threshold
evenlim = int(dict['def']) # numero limite de meses para contar un evento de sequia
mask = int(dict['lsm']) #define si se aplica mascara de tierra (1) o no (0)
extdrou = float(dict['extdrou'])
sevdrou = float(dict['sevdro'])
moddrou = float(dict['moddro'])


#Por el momento desabilitar esta linea y coger siempre 1
#set etp = $argv[20] # 0=etp_original->03312, 1=etp_thornthwaite->03312T
etp = 1


root_dir = os.getcwd()

home = root_dir+'/child-processes/data'
# home = 'child-processes/data'
workdir = '/'.join([home,'data4drought'])
datadir = 'data'
progdir = 'programs'
utildir = 'utils'
outdir = 'outputs'

#fijando irmon,ermon, iamon,eamon
irmon = 1
ermon = 12
iamon = 1
eamon = 12


# se define el valor missing que tendran todos los files que se procesen
undef = -999999

# if dindex==1: #spi
#     ifilename = 'pm.05216'
# else: 
#     ifilename = 'pm.63312'


# ifile = '/'.join([workdir, datadir, modeldir2, ifilename])


# # SPI NC_FILE
# spi_nc_name = 'pm.05216'
spi_nc_name = 'pr-ACCESS-ESM1'

# # SPEI NC_FILE
spei_nc_name = 'pm.063312_wtb'


if dindex==1: #spi
    if mask==0:
        ifilename = spi_nc_name
    else:
        ifilename = spi_nc_name+'_mask'
else: 
    if mask==0:
        ifilename = spei_nc_name
    else:
        ifilename = spei_nc_name+'_mask'


ifile = '/'.join([workdir, datadir, modeldir2, ifilename])

if not os.path.exists(ifile+'.nc'):
    subprocess.run('cat {} >> {}'.format('/'.join([workdir, datadir, modeldir2, 'nc_file.parta*']), ifile+'.nc.tar.gz'), shell=True, capture_output=True)
    subprocess.run('tar -xzvf {} -C {}'.format(ifile+'.nc.tar.gz','child-processes/data/data4drought/data/echam5/'), shell=True, capture_output=True)
    subprocess.run('rm {}'.format(ifile+'.nc.tar.gz'), shell=True, capture_output=True)


# targzfiles = glob.glob('/'.join([workdir, datadir, modeldir2,'*.tar.gz']))
# for file in targzfiles:
    # subprocess.run('tar -xzvf {} -C {}'.format(file,'child-processes/data/data4drought/data/echam5/'), shell=True, capture_output=True)


# namevar_cdo = subprocess.run('cdo -s showname {}.nc'.format(ifile), shell=True, encoding='utf-8', capture_output=True).stdout
# # ##print('ifile',ifile)
# # ##print('namevar_cdo',namevar_cdo)
# namevar = namevar_cdo.split(' ')[1].split('\n')[0]

# #print('aqui',ifile, namevar)
# exit()

# def find_files(filename, search_path):
#    result = []

# # Wlaking top-down from the root
#    for root, dir, files in os.walk(search_path):
#       if filename in files:
#          result.append(os.path.join(root, filename))
#    return result

# result = find_files(".bashrc","/home")

# #print(result)

# #print(find_files("g++","/usr"))
# exit()

# #print('********* nc-config *********')
# #print()
# subprocess.run('nc-config --all', shell=True, capture_output=True)
# #print()
# #print('********* nc-config *********')


def install_cdo():
    os.chdir('/'.join([root_dir,'child-processes/cdo']))
    # subprocess.run('rm cdo.log', shell=True)
    sh_file = 'install.sh'
    subprocess.run('chmod a+x {}'.format(sh_file), shell=True)
    subprocess.run(f'./{sh_file} >> cdo.log', shell=True)
    # #print('./configure success!')
    # #print("CDO was instaled success!")

# install_cdo()

# dir_child_process = '/'.join([home_server, 'child-processes'])
# cdo_dir = '/'.join([dir_child_process,'cdo/bin/cdo'])
# # # local_cdo_dir = '/home/adrianfb/cdo_install/cdo-1.9.1/local/bin/cdo'
# # # #print(cdo_dir)
# # #print()
cdo_version = subprocess.run('cdo --version', shell=True, encoding='utf-8', capture_output=True)
# # cdo_version = subprocess.run('{} --version'.format(local_cdo_dir), shell=True, encoding='utf-8', capture_output=True)
#print('line00', cdo_version.stdout.split('\n')[0])
#print('line11',cdo_version.stdout.split('\n')[1])

exit()

os.chdir(home_server)

subprocess.run('git config http.postBuffer 1048576000', shell=True)
subprocess.run('git config user.email "adrianfuentesbarrios@gmail.com"', shell=True)
subprocess.run('git config user.name "adrianfb94"', shell=True)


# #print('estoy haciendo el init')
# subprocess.run('git init', shell=True, capture_output=False)
# # os.system('git init')
# #print('ya hice el init')

# #print("estoy haciendo el add zlib")
# subprocess.run('git add -A {}'.format('/'.join(['child-processes','zlib-1.2.8'])), shell=True, capture_output=False)
# # os.system('git add file_added_from_server.txt')
# #print('ya hice el add')

# #print("estoy haciendo el commit zlib")
# subprocess.run('git commit -m "added zlib from SERVER"', shell=True, capture_output=False)
# # os.system('git commit -m "added from SERVER"')
# #print('ya hice el commit')


# #print("estoy haciendo el add hdf5")
# subprocess.run('git add -A {}'.format('/'.join(['child-processes','hdf5-1.8.13'])), shell=True, capture_output=False)
# # os.system('git add file_added_from_server.txt')
# #print('ya hice el add')

# #print("estoy haciendo el commit hdf5")
# subprocess.run('git commit -m "added hdf5 from SERVER"', shell=True, capture_output=False)
# # os.system('git commit -m "added from SERVER"')
# #print('ya hice el commit')


# #print("estoy haciendo el add netcdf")
# subprocess.run('git add -A {}'.format('/'.join(['child-processes','netcdf-c-4.5.0'])), shell=True, capture_output=False)
# # os.system('git add file_added_from_server.txt')
# #print('ya hice el add')

# #print("estoy haciendo el commit netcdf")
# subprocess.run('git commit -m "added netcdf from SERVER"', shell=True, capture_output=False)
# # os.system('git commit -m "added from SERVER"')
# #print('ya hice el commit')



# #print("estoy haciendo el add cdo")
# subprocess.run('git add -A {}'.format('/'.join(['child-processes','cdo-1.9.1'])), shell=True, capture_output=True)
# # os.system('git add file_added_from_server.txt')
# #print('ya hice el add')

# #print("estoy haciendo el commit cdo")
# subprocess.run('git commit -m "added cdo from SERVER"', shell=True, capture_output=True)
# # os.system('git commit -m "added from SERVER"')
# #print('ya hice el commit')


# #print("estoy haciendo el add cdo_installed")
# subprocess.run('git add -A {}'.format('/'.join(['child-processes','cdo_installed'])), shell=True, capture_output=False)
# # os.system('git add file_added_from_server.txt')
# #print('ya hice el add')

# #print("estoy haciendo el commit cdo_installed")
# subprocess.run('git commit -m "added cdo_installed from SERVER"', shell=True, capture_output=False)
# # os.system('git commit -m "added from SERVER"')
# #print('ya hice el commit')


#print("estoy haciendo el add -A")
subprocess.run('git add -A', shell=True, capture_output=False)
# os.system('git add file_added_from_server.txt')
#print('ya hice el add -A')

#print("estoy haciendo el commit cdo_installed")
subprocess.run('git commit -m "added from SERVER"', shell=True, capture_output=False)
# os.system('git commit -m "added from SERVER"')
#print('ya hice el commit')


# #print("estoy haciendo el pull")
# subprocess.run('git pull -v https://{}@github.com/{}/{}.git HEAD:main'.format(token, username, repo_name), shell=True, capture_output=False)
# # os.system(f'git push https://{token}@github.com/{username}/{repo_name}.git HEAD:main')
# #print('ya hice el pull')



#print("estoy haciendo el push")
subprocess.run(f'git push https://{token}@github.com/{username}/{repo_name}.git HEAD:main', shell=True, capture_output=False)
# subprocess.run(f'git push https://github.com/{username}/{repo_name}.git HEAD:main', shell=True, capture_output=False)
#print('ya hice el push')



exit()
cdo_dir = '/'.join([root_dir,'child-processes/cdo-1.9.1/local/bin/cdo'])
# local_cdo_dir = '/home/adrianfb/cdo_install/cdo-1.9.1/local/bin/cdo'

# #print(cdo_dir)

cdo_version = subprocess.run('{} --version'.format(cdo_dir), shell=True, encoding='utf-8', capture_output=True)
# cdo_version = subprocess.run('{} --version'.format(local_cdo_dir), shell=True, encoding='utf-8', capture_output=True)
#print(cdo_version.stderr.split('\n')[0])
#print(cdo_version.stderr.split('\n')[1])


# from cdo import *
# cdo = Cdo()
# #print('\n')
# #print('cdo.version(): {}'.format(cdo.version()))
# #print()
# nvar = cdo.showname(input=ifile+'.nc')
# #print("namevar: {}".format(nvar))
# #print()
# grid = cdo.griddes(input=ifile+'.nc')
# #print("griddes: {}".format(grid))
# #print()
# ntime = cdo.ntime(input=ifile+'.nc')
# #print("ntime: {}".format(ntime))

exit()

# tempvar = subprocess.run('cdo -s showname {}.nc'.format(ifile), shell=True, encoding='utf-8', capture_output=True).stdout
namevar = cdo.showname(input=ifile)
#print('namevar: {}'.format(namevar))
exit()


##print('cdo ifile', ifile)
grid = subprocess.run('cdo -s griddes {}.nc'.format(ifile), shell=True, encoding='utf-8', capture_output=True).stdout

#0 #
#1 # gridID 1
#2 #
#3 gridtype  = lonlat
#4 gridsize  = 28085
#5 xsize     = 205
#6 ysize     = 137
#7 xname     = lon
#8 xlongname = "Longitude"
#9 xunits    = "degrees_east"
#10 yname     = lat
#11 ylongname = "Latitude"
#12 yunits    = "degrees_north"
#13 xfirst    = -94.9
#14 xinc      = 0.22
#15 yfirst    = 0.0199998915
#16 yinc      = 0.220001
#17 scanningMode = 64



gridlines = grid.split('\n')
#print(gridlines)
exit()

nx = int(gridlines[5].split('=')[1])
# ##print('nx = {}'.format(nx))

ny = int(gridlines[6].split('=')[1])
# ##print('ny = {}'.format(ny))

xfirst = float(gridlines[13].split('=')[1])
# ##print('xfirst = {}'.format(xfirst))

xinc = float(gridlines[14].split('=')[1])
# ##print('xinc = {}'.format(xinc))

yfirst = float(gridlines[15].split('=')[1])
# ##print('yfirst = {}'.format(yfirst))

yinc = float(gridlines[16].split('=')[1])
# ##print('yinc = {}'.format(yinc))


ntime = int(subprocess.run('cdo -s ntime {}.nc'.format(ifile), shell=True, encoding='utf-8', capture_output=True).stdout)
# ##print('ntime = {}'.format(ntime))


args = ['INIT_0', dindex, modeldir, workdir, datadir, nx, xfirst, xinc, ny, 
        yfirst, yinc, ntime, undef, ilon, elon, ilat, elat, avg, 
        iryear, eryear, irmon, ermon, iayear, eayear, iamon, eamon, 
        spilen, spilim, evenlim, etp, mask, extdrou, sevdrou, moddrou
        ]



# gradscmdline = [grads, mode, file]

# for arg in args:
#     gradscmdline.append(arg)

# ##print()
# for i,arg in enumerate(args):
#     ##print(i,arg, type(arg))


'''
Aqui comienza GRADS
'''
##print()

def namemonth(iimes,a):
    nmes = []
    nmes.append('jan JAN Enero')
    nmes.append('feb FEB Febrero')
    nmes.append('mar MAR Marzo')
    nmes.append('apr APR Abril')
    nmes.append('may MAY Mayo')
    nmes.append('jun JUN Junio')
    nmes.append('jul JUL Julio')
    nmes.append('aug AUG Agosto')
    nmes.append('sep SEP Septiembre')
    nmes.append('oct OCT Octubre')
    nmes.append('nov NOV Noviembre')
    nmes.append('dec DEC Diciembre')
    namemes=nmes[iimes-1].split(' ')[a-1]
    return(namemes)




didx=args[1]
model=args[2] #subdirectorio del experimento
modeldir=models[model]# the getmodel function is at the end of this script
# if(didx==1): 
#     dindex='spi'
#     var1='precip'

# else:
#     dindex='spei'
#     var1='wtb'

if(didx==1): 
    dindex='spi'
    var1='precip'

else:
    dindex='spei'
    var1='wtb'






_workdir=args[3] #directorio general de trabajo "home/<username>/<directorio>"
datadir=args[4] #nombre del directorio donde se ubica <modeldir>
nx=args[5] #numero de puntos por las xs (longitudes)
xfirst=args[6] #x de inicio
xinc=args[7] #incremento por las x
ny=args[8] #numero de puntos por las ys (latitudes)
yfirst=args[9] #y de inicio
yinc=args[10] #incremento por las y (si viene negativo, se hace positivo)
if(yinc<0):
    yinc=yinc*(-1)
else:
    yinc=yinc
time=args[11] #numero de tiempos (se extrae del CDO)
undef=args[12] #valor missing del fichero
ilon=args[13] #longitud inicial
elon=args[14] #longitud final
ilat=args[15] #latitud inicial
elat=args[16] #latitud final
avg=args[17] #0->all_grid_points,1->area_average,2->latitude_average,4->longitude_average
iryear=args[18] #anno de inicio del periodo de referencia
eryear=args[19] #anno final del periodo de referencia
irmon=args[20] #mes de inicio (creo que siempre sera 1) periodo de referencia
ermon=args[21] #mes final (creo que siempre sera 12) periodo de referencia
iayear=args[22] # anno de inicio del periodo de analisis
eayear=args[23] #anno final de periodo de analisis
iamon=args[24] #mes de inicio periodo analisis
eamon=args[25] #mes final periodo de analisis
spilen=args[26] #cantidad de meses para spi o spei
spilim=args[27] #valor limite de spi o spei para la contabilidad
evenlim=args[28] #valor del numero de meses para contar periodos de sequia
pet=args[29] # define que etp usa. 0->lee file pm.03312.nc, 1->lee file pm.03312T.nc 
mask=args[30] # define si se usa mascara de tierra (1) o no se usa (0)
extdrou=args[31] #extreme drought (spi or spei) value
sevdrou=args[32] #severe drought (spi or spei) value
moddrou=args[33] #moderate drought (spi or spei) value

avgo=avg

lsm=mask
##print(modeldir,ilon,ilat,elon,elat,iryear,eryear,iayear,eayear,dindex+str(spilen),dindex,spilim,dindex,evenlim,extdrou,sevdrou,moddrou,var1)
##print(dindex+str(spilen))



if(modeldir=='crudata' or modeldir=='udeldata'or modeldir=='gpccdata'):
    lsm=0

if(modeldir=='cmapdata'):
    lsm=1

irmon=namemonth(irmon,2)
ermon=namemonth(ermon,2)
iamon=namemonth(iamon,2)
eamon=namemonth(eamon,2)


irdate=irmon+str(iryear)
erdate=ermon+str(eryear)
iadate=iamon+str(iayear)
eadate=eamon+str(eayear)

# if(dindex=='spi'):
#     if(lsm==0):
#         filedir = ('/'.join([_workdir, datadir, modeldir, 'pm.05216.nc']))
#     else:
#         filedir = ('/'.join([_workdir, datadir, modeldir, 'pm.05216_mask.nc']))


#     ofile1='_'.join([modeldir, str(ilon), str(ilat), str(elon), str(elat), str(iryear), 
#                 str(eryear), str(iayear), str(eayear),dindex+str(spilen),dindex,
#                 str(spilim), dindex, str(evenlim), str(extdrou), str(sevdrou),
#                 str(moddrou), var1])
#     voper=var1
# else:
#     if(lsm==0):
#         filedir = ('/'.join([_workdir, datadir, modeldir, 'pm.63312.nc']))
#     else:
#         filedir = ('/'.join([_workdir, datadir, modeldir, 'pm.63312_mask.nc']))

#     ofile1='_'.join([modeldir, str(ilon), str(ilat), str(elon), str(elat), str(iryear), 
#                 str(eryear), str(iayear), str(eayear),dindex+str(spilen),dindex,
#                 str(spilim), dindex, str(evenlim), str(extdrou), str(sevdrou),
#                 str(moddrou), var1])

#     voper=var1

if(dindex=='spi'):
    if(lsm==0):
        filedir = ('/'.join([_workdir, datadir, modeldir, spi_nc_name+'.nc']))
    else:
        filedir = ('/'.join([_workdir, datadir, modeldir, spi_nc_name+'_mask.nc']))


    ofile1='_'.join([modeldir, str(ilon), str(ilat), str(elon), str(elat), str(iryear), 
                str(eryear), str(iayear), str(eayear),dindex+str(spilen),dindex,
                str(spilim), dindex, str(evenlim), str(extdrou), str(sevdrou),
                str(moddrou), var1])
    voper=var1
else:
    if(lsm==0):
        filedir = ('/'.join([_workdir, datadir, modeldir, spei_nc_name+'.nc']))
    else:
        filedir = ('/'.join([_workdir, datadir, modeldir, spei_nc_name+'_mask.nc']))

    ofile1='_'.join([modeldir, str(ilon), str(ilat), str(elon), str(elat), str(iryear), 
                str(eryear), str(iayear), str(eayear),dindex+str(spilen),dindex,
                str(spilim), dindex, str(evenlim), str(extdrou), str(sevdrou),
                str(moddrou), var1])

    voper=var1




def closest(lst, K):
     
     lst = np.asarray(lst)
     idx = (np.abs(lst - K)).argmin()
     return lst[idx], idx+1






# def save_mean(avg, minlat, maxlat, minlon, maxlon, ini, end, axis=None):
#     list_variable = []

#     if avg != 0:
#         if avg==1:

#             ##print('AVG = {}'.format(avg))
#             ##print('ini {} end {}'.format(ini, end))
#             ##print('minlon {} maxlon {}'.format(minlon, maxlon))
#             ##print('ini {} end {}'.format(0, nc_variable.shape[0]))


#             lon, lat = np.meshgrid(nc_lon, nc_lat)

#             Lat_cond = (lat >= minlat) & (lat <= maxlat)
#             Lon_cond = (lon >= minlon) & (lon <= maxlon)

#             lat_lon = Lat_cond & Lon_cond

#             num_lat = np.where(Lat_cond[:,0] == True)
#             num_lon = np.where(Lon_cond[0,:] == True)

#             dim_lat, dim_lon = len(num_lat[0]), len(num_lon[0])

#             lat = lat[lat_lon].reshape((dim_lat, dim_lon))
#             lon = lon[lat_lon].reshape((dim_lat, dim_lon))

#             mask_nc_variable = nc_variable==undef
#             nc_variable_masked = ma.masked_array(nc_variable, mask=mask_nc_variable)


#             # list_variable = []
#             if axis==None:
#                 for k in range(ini-1,end):
#                     variable = nc_variable_masked[k,:,:]
#                     variable = variable[lat_lon].reshape((dim_lat, dim_lon))
#                     list_variable.append(variable.mean())                    
#                 list_variable = np.asarray(list_variable)
#                 # np.savetxt('./data4drought/outputs/prueba/local/127.0.0.1/output_python.txt', list_variable)
#                 # #print('file saved')
#                 # #print(list_variable.shape)
#                 # return list_variable
#             elif axis==0:
#                 # ##print('var not found')
#                 for k in range(ini-1,end):
#                 # for k in range(0,nc_variable.shape[0]):
#                     variable = nc_variable_masked[k,:,:]
#                     variable = variable[lat_lon].reshape((dim_lat, dim_lon))
#                     list_variable.append(variable.mean(axis=0).mean())
#                 list_variable = np.asarray(list_variable)
#                 # np.savetxt('./data4drought/outputs/prueba/output_python.txt', list_variable)
#                 # ##print('file saved with lon')
#             elif axis==1:
#                 for k in range(ini-1,end):
#                 # for k in range(0,nc_variable.shape[0]):

#                     variable = nc_variable_masked[k,:,:]
#                     variable = variable[lat_lon].reshape((dim_lat, dim_lon))
#                     list_variable.append(variable.mean(axis=1).mean())
#                 list_variable = np.asarray(list_variable)
#                 # np.savetxt('./data4drought/outputs/prueba/output_python.txt', list_variable)
#                 # ##print('file saved with lat')

#             return list_variable
#         else:
#             #print('No code for avg={}'.format(avg))
#             exit()

#     else:
#         ##print('AVG = 0'.format(avg))
#         ##print('ini {} end {}'.format(ini, end))
#         ##print('minlon {} maxlon {}'.format(minlon, maxlon))
#         ##print('ini {} end {}'.format(ini-1, end))

#         lon, lat = np.meshgrid(nc_lon, nc_lat)

#         Lat_cond = (lat >= minlat) & (lat <= maxlat)
#         Lon_cond = (lon >= minlon) & (lon <= maxlon)

#         lat_lon = Lat_cond & Lon_cond

#         num_lat = np.where(Lat_cond[:,0] == True)
#         num_lon = np.where(Lon_cond[0,:] == True)

#         dim_lat, dim_lon = len(num_lat[0]), len(num_lon[0])

#         lat = lat[lat_lon].reshape((dim_lat, dim_lon))
#         lon = lon[lat_lon].reshape((dim_lat, dim_lon))
#         for k in range(ini-1,end):
#         # for k in range(0,nc_variable.shape[0]):
#             variable = nc_variable[k,:,:]
#             variable = variable[lat_lon].reshape((dim_lat, dim_lon))
#             list_variable.append(variable)
#         list_variable = np.asarray(list_variable)
#         ##print('shape list_variable', list_variable.shape)
#         ##print('min list_variable', list_variable.min())
#         ##print('max list_variable', list_variable.max())
#         # np.savetxt('./data4drought/outputs/prueba/output_python.txt', list_variable)
#         # ##print('file saved')
#         return list_variable

def save_mean(avg, minlat, maxlat, minlon, maxlon, ini, end, axis=None):
    weighted = True

    list_variable = []

    if avg != 0:
        if avg==1:

            mask_nc_variable = nc_variable==undef
            nc_variable_masked = ma.masked_array(nc_variable, mask=mask_nc_variable)

            lon, lat = np.meshgrid(nc_lon, nc_lat)

            Lat_cond = (lat >= minlat) & (lat <= maxlat)
            Lon_cond = (lon >= minlon) & (lon <= maxlon)

            lat_lon = Lat_cond & Lon_cond

            num_lat = np.where(Lat_cond[:,0] == True)
            num_lon = np.where(Lon_cond[0,:] == True)

            dim_lat, dim_lon = len(num_lat[0]), len(num_lon[0])

            lat = lat[lat_lon].reshape((dim_lat, dim_lon))
            lon = lon[lat_lon].reshape((dim_lat, dim_lon))

            if axis==None:
                if not weighted:
                    for k in range(ini-1,end):
                        variable = nc_variable_masked[k,:,:]
                        variable = variable[lat_lon].reshape((dim_lat, dim_lon))
                        list_variable.append(variable.mean())                    

                    list_variable = np.asarray(list_variable)
                    # np.savetxt('./data4drought/outputs/prueba/local/127.0.0.1/output_python_no_weighted.txt', list_variable)

                else:


                    # ds = xr.load_dataset(filedir)
                    # if (ds.lon.values.min()>0 and ds.lon.values.max()>0):
                    #     ds['lon'] = ds['lon']-360
                    # mask_lon = (ds.lon >= minlon) & (ds.lon <= maxlon)
                    # mask_lat = (ds.lat >= minlat) & (ds.lat <= maxlat)   
                    # cropped_ds = ds.where(mask_lon & mask_lat, drop=True)

                    mask_lon = (nc_lon >= minlon) & (nc_lon <= maxlon)
                    mask_lat = (nc_lat >= minlat) & (nc_lat <= maxlat)   
                    cropped_ds = nc_all.where(mask_lon & mask_lat, drop=True)

                    weights = np.cos(np.deg2rad(cropped_ds.lat))
                    weights.name = "weights"

                    data = cropped_ds[var1]
                    data_weighted = data.weighted(weights)
                    weighted_mean = data_weighted.mean(("lon", "lat"))
                    list_variable = np.asarray(weighted_mean[ini-1:end].values)
                    # np.savetxt('./data4drought/outputs/prueba/local/127.0.0.1/output_python_weighted.txt', list_variable)



            # elif axis==0:
            #     for k in range(ini-1,end):
            #         variable = nc_variable_masked[k,:,:]
            #         variable = variable[lat_lon].reshape((dim_lat, dim_lon))
            #         list_variable.append(variable.mean(axis=0).mean())

            #     list_variable = np.asarray(list_variable)

            # elif axis==1:
            #     for k in range(ini-1,end):
            #         variable = nc_variable_masked[k,:,:]
            #         variable = variable[lat_lon].reshape((dim_lat, dim_lon))
            #         list_variable.append(variable.mean(axis=1).mean())

            #     list_variable = np.asarray(list_variable)

            return list_variable
        else:
            # #print('no tengo codigo para avg = {}'.format(avg))
            exit()
    
    else:
    
        lon, lat = np.meshgrid(nc_lon, nc_lat)
        Lat_cond = (lat >= minlat) & (lat <= maxlat)
        Lon_cond = (lon >= minlon) & (lon <= maxlon)

        lat_lon = Lat_cond & Lon_cond

        num_lat = np.where(Lat_cond[:,0] == True)
        num_lon = np.where(Lon_cond[0,:] == True)

        dim_lat, dim_lon = len(num_lat[0]), len(num_lon[0])

        lat = lat[lat_lon].reshape((dim_lat, dim_lon))
        lon = lon[lat_lon].reshape((dim_lat, dim_lon))

        for k in range(ini-1,end):
            variable = nc_variable[k,:,:].values
            variable = variable[lat_lon].reshape((dim_lat, dim_lon))
            list_variable.append(variable)

        list_variable = np.asarray(list_variable)

        return list_variable







# def read_from_grb(file):
#     f = open(file, 'rb')
#     data = f.read()
#     f.close()
#     assert len(data)%4==0
#     count = len(data)//4
#     result = struct.unpack('<{0}f'.format(count), data)
#     return result



def righcoor(avg,ilat,elat,ilon,elon):

    if(avg!=0 and ilat==elat and ilon==elon):
        avg=0
    if(avg==1 and ilat==elat and ilon!=elon):
      avg=3
    if(avg==1 and ilat!=elat and ilon==elon):
      avg=2
    if(avg==2 and ilat==elat and ilon!=elon):
      avg=3
    if(avg==3 and ilat!=elat and ilon==elon):
      avg=2

    if(avg==0):
        yfirst, y1 = closest(nc_lat, ilat)
        ylast, y2 = closest(nc_lat, elat)
        ny=y2-y1+1

        xfirst, x1 = closest(nc_lon, ilon)
        xlast, x2 = closest(nc_lon, elon)
        nx=x2-x1+1

    if(avg==1):
        xfirst=(ilon+elon)/2
        yfirst=(ilat+elat)/2
        xlast=xfirst
        ylast=yfirst
        nx=1
        ny=1

    if(avg==2):
        xfirst, x1 = closest(nc_lon, ilon)
        xlast, x2 = closest(nc_lon, elon)
        nx=x2-x1+1
        yfirst=(ilat+elat)/2
        ylast=yfirst
        ny=1
    
    if(avg==3):
        yfirst, y1 = closest(nc_lat, ilat)
        ylast, y2 = closest(nc_lat, elat)
        ny=y2-y1+1
        xfirst=(ilon+elon)/2
        xlast=xfirst
        nx=1

    return xfirst,yfirst,nx,ny,xlast,ylast,avg

# def searchtime_date(time_ref):
#     # dates = np.array(nc_time)
#     dates = nc_time_array
#     # dates = num2date(dates, units=nc_time.units)
#     dates = num2date(dates, units=nc_time_units)
#     for i,date in enumerate(dates):
#         date_tmp = str(date)
#         year = date_tmp.split(' ')[0].split('-')[0]
#         month = date_tmp.split(' ')[0].split('-')[1]
#         current_date = namemonth(int(month),2)+year
#         if current_date == time_ref:
#             K = i+1
#     return K

def searchtime_date(time_ref):

    dates = nc_time_array
    # dates = num2date(dates, units=nc_time_units)

    for i,date in enumerate(dates):
        date_tmp = str(date)
        year = date_tmp.split(' ')[0].split('-')[0]
        month = date_tmp.split(' ')[0].split('-')[1]

        current_date = namemonth(int(month),2)+year
        if current_date == time_ref:
            K = i+1
    return K



# def searchtime_pos(pos_i, pos_f):
#     # dates = np.array(nc_time)
#     dates = nc_time_array
#     # dates = num2date(dates, units=nc_time.units)
#     dates = num2date(dates, units=nc_time_units)
#     date_i = str(dates[pos_i])
#     year = date_i.split(' ')[0].split('-')[0]
#     month = date_i.split(' ')[0].split('-')[1]
#     day = date_i.split(' ')[0].split('-')[2]
#     hour = date_i.split(' ')[1].split(':')[0]+'Z'
#     date_i = hour+day+namemonth(int(month),2)+year
#     date_f = str(dates[pos_f])
#     year = date_f.split(' ')[0].split('-')[0]
#     month = date_f.split(' ')[0].split('-')[1]
#     day = date_f.split(' ')[0].split('-')[2]
#     hour = date_f.split(' ')[1].split(':')[0]+'Z'
#     date_f = hour+day+namemonth(int(month),2)+year

#     # ##print('Time = {} to {} T = {} to {}'.format(date_i, date_f, pos_i+1, pos_f+1))

#     time = pos_f - pos_i + 1


#     return date_i, date_f, time

def searchtime_pos(pos_i, pos_f):
    dates = nc_time_array
    date_i = str(dates[pos_i])
    year = date_i.split(':')[0].split('-')[0]
    month = date_i.split(':')[0].split('-')[1]
    day = date_i.split(':')[0].split('-')[2][0:2]
    hour = date_i.split(':')[0].split('-')[2][3:5]+'Z'
    date_i = hour+day+namemonth(int(month),2)+year

    date_f = str(dates[pos_f])
    year = date_f.split(':')[0].split('-')[0]
    month = date_f.split(':')[0].split('-')[1]
    day = date_f.split(':')[0].split('-')[2][0:2]
    hour = date_f.split(':')[0].split('-')[2][3:5]+'Z'
    date_f = hour+day+namemonth(int(month),2)+year

    time = pos_f - pos_i + 1

    return date_i, date_f, time



def tdates(t_in):
    t_out = []
    i=0
    while(i<4):
        t_out.append(searchtime_date(t_in[i]))
        i+=1
    if(t_out[0]>t_out[2]):
        tmin=t_out[2]
    else:
        tmin=t_out[0]
    if(t_out[1]<t_out[3]):
        tmax=t_out[3]
    else:
        tmax=t_out[1]
    return tmin, tmax, t_out


def avg_name(avg):
    if(avg==1):
        average='average'
    if(avg==0):
        average='gridboxes'
    if(avg==2):
        average='longcross' # se promedian las latitudes
    if(avg==3):
        average='latcross' # se promedian las longitudes
    
    return(average)



# def operation(tmin, tmax, namefile):
#     ##print('operation avg es ', avg)
#     ##print('operation namefile es ', namefile)
#     ##print('operation tmin es ', tmin)
#     ##print('operation tmax es ', tmax)
#     '''
#     ****************************************************************
#     * Return the command line that will be executed to produce area average
#     * gridbox, lat average or long average.
#     *
#     * This dynamical script function is called by spi_spei.gs GrADS script
#     * 
#     * Writen by: Abel Centella
#     *            October, 2014
#     ****************************************************************
#     '''
#     option = ''
#     if(avg==1):
#         var_array = save_mean(avg, ini=tmin, end=tmax, minlat=ilat, minlon=ilon, maxlat=elat, maxlon=elon, axis=None)
#         var_bytes = bytes()
#         for val in var_array:
#             var_bytes += struct.pack('f', val)
#         file_grb = open(namefile+'.grb', 'wb')
#         file_grb.write(var_bytes)
#         file_grb.close()
#         ##print('avg={} => grb saved'.format(avg))
#         if(option==''):
#             option='Area_Average'
#     if(avg==0):
#         var_array = save_mean(avg, ini=tmin, end=tmax, minlat=yfirst, minlon=xfirst, maxlat=ylast, maxlon=xlast, axis=None)
#         ##print(10*'@')
#         ##print(var_array.shape)
#         ##print(type(var_array))
#         ##print(10*'@')

#         # var_bytes = bytes()
#         # for val in var_array:
#         #     ##print('{}/{}'.format(list(var_array).index(val), len(var_array)))
#         #     var_bytes += struct.pack('f', val)
#         file_grb = open(namefile+'.grb', 'wb')
#         file_grb.write(var_array.astype('f').tobytes())
#         # file_grb.write(struct.pack('{}f'.format([var for var in var_array]), var_array))
#         file_grb.close()
#         ##print('acabo de guardar grb para avg=0')

#         if(option=='' and ilon!=elon and ilat!=elat):
#             option='Grid_Point_Values'
#         if(option=='' and ilon==elon and ilat!=elat):
#             option='Longitude_Cross_Section'
#         if(option=='' and ilon!=elon and ilat==elat):
#             option='Latitude_Cross_Section'
#         if(option=='' and ilon==elon and ilat==elat):
#             option='Grid_Point_Value'
#     if(avg==2):
#         # operation='const(ave('var',lat='ilat',lat='elat'),'undef',-u)'
#         var_array = save_mean(avg, ini=tmin, end=tmax, minlat=ilat, minlon=ilon, maxlat=elat, maxlon=elon, axis=1)
#         var_bytes = bytes()
#         for val in var_array:
#             var_bytes += struct.pack('f', val)
#         file_grb = open(namefile+'.grb', 'wb')
#         file_grb.write(var_bytes)
#         file_grb.close()

#         ##print('avg={} => grb saved'.format(avg))
#         # 'set lon 'xfirst' 'xlast
#         # 'set lat 'yfirst' 'ylast
#         if(option==''):
#             option='Latitudinal_Average'
#     if(avg==3):
#         # operation='const(ave('var',lon='ilon',lon='elon'),'undef',-u)'
#         var_array = save_mean(avg, ini=tmin, end=tmax, minlat=ilat, minlon=ilon, maxlat=elat, maxlon=elon, axis=0)

#         var_bytes = bytes()
#         for val in var_array:
#             var_bytes += struct.pack('f', val)
#         file_grb = open(namefile+'.grb', 'wb')
#         file_grb.write(var_bytes)
#         file_grb.close()



#         ##print('avg={} => grb saved'.format(avg))
#         # 'set lon 'xfirst' 'xlast
#         # 'set lat 'yfirst' 'ylast
#         if(option==''):
#             option='Longitudinal_Average'

#     average=avg_name(avg)

#     return option, average

def operation(tmin, tmax, namefile):

    option = ''
    if(avg==1):

        var_array = save_mean(avg, ini=tmin, end=tmax, minlat=ilat, minlon=ilon, maxlat=elat, maxlon=elon, axis=None)

        var_bytes = bytes()
        for val in var_array:
            var_bytes += struct.pack('f', val)

        file_grb = open(namefile+'.grb', 'wb')
        file_grb.write(var_bytes)
        file_grb.close()

        if(option==''):
            option='Area_Average'

    if(avg==0):

        var_array = save_mean(avg, ini=tmin, end=tmax, minlat=yfirst, minlon=xfirst, maxlat=ylast, maxlon=xlast, axis=None)

        file_grb = open(namefile+'.grb', 'wb')
        file_grb.write(var_array.astype('f').tobytes())
        file_grb.close()

        if(option=='' and ilon!=elon and ilat!=elat):
            option='Grid_Point_Values'
        if(option=='' and ilon==elon and ilat!=elat):
            option='Longitude_Cross_Section'
        if(option=='' and ilon!=elon and ilat==elat):
            option='Latitude_Cross_Section'
        if(option=='' and ilon==elon and ilat==elat):
            option='Grid_Point_Value'

    if(avg==2):

        var_array = save_mean(avg, ini=tmin, end=tmax, minlat=ilat, minlon=ilon, maxlat=elat, maxlon=elon, axis=1)

        var_bytes = bytes()
        for val in var_array:
            var_bytes += struct.pack('f', val)

        file_grb = open(namefile+'.grb', 'wb')
        file_grb.write(var_bytes)
        file_grb.close()

        if(option==''):
            option='Latitudinal_Average'

    if(avg==3):

        var_array = save_mean(avg, ini=tmin, end=tmax, minlat=ilat, minlon=ilon, maxlat=elat, maxlon=elon, axis=0)

        var_bytes = bytes()
        for val in var_array:
            var_bytes += struct.pack('f', val)

        file_grb = open(namefile+'.grb', 'wb')
        file_grb.write(var_bytes)
        file_grb.close()

        if(option==''):
            option='Longitudinal_Average'

    average=avg_name(avg)

    return option, average





def write_ctl(path, fname,nx,xfirst,xinc,ny,yfirst,yinc,time,date,undef,var,indx):
    ##print('this is write_ctl function', fname)
    # ##print('xdef {} linear {} {}\n'.format(nx, xfirst, xinc))
    fnamectl = fname + '.ctl'
    octl = '/'.join([path, fnamectl])
    # ##print(path, fname, octl)
    with open(octl, 'w') as file_ctl:
        file_ctl.write('Dset ^{}.grb\n'.format(fname))
        file_ctl.write('title drought index data input\n')
        file_ctl.write('undef {}\n'.format(undef))
        file_ctl.write('xdef {} linear {} {}\n'.format(nx, xfirst, xinc))
        file_ctl.write('ydef {} linear {} {}\n'.format(ny, yfirst, yinc))
        file_ctl.write('zdef 1 linear 1 1\n')
        file_ctl.write('tdef {} linear {} 1mo\n'.format(time, date))
        file_ctl.write('vars 1\n')
        file_ctl.write('{}{} 0 {} **\n'.format(var, indx, undef))
        file_ctl.write('endvars\n')
        file_ctl.write(' \n')
        file_ctl.close()
    
    ##print('ctl file [{}] was created'.format(fnamectl))
    ##print()

    return



# nc_lat, nc_lon, nc_variable, nc_time_array, nc_time_units = get_nc(var1, filedir ,undef)
# xfirst,yfirst,nx,ny,xlast,ylast,avg = righcoor(avg,ilat,elat,ilon,elon)
# ##print('xfirst:{}, yfirst:{}, nx:{}, ny:{}, xlast:{}, ylast:{}, avg:{}'.format(xfirst,yfirst,nx,ny,xlast,ylast,avg))
# ##print('irdate={}, erdate={}, iadate={}, eadate={}'.format(irdate, erdate, iadate, eadate))
# t_in = [irdate, erdate, iadate, eadate]
# tmin, tmax, [t1, t2, t3, t4] = tdates(t_in)


# nc_lat, nc_lon, nc_variable, nc_time_array, nc_time_units = get_nc(var1, filedir, undef)
nc_all, nc_lat, nc_lon, nc_variable, nc_time_array = get_nc(var1, filedir, undef)


xfirst,yfirst,nx,ny,xlast,ylast,avg = righcoor(avg,ilat,elat,ilon,elon)
t_in = [irdate, erdate, iadate, eadate]
tmin, tmax, [t1, t2, t3, t4] = tdates(t_in)


itr = t1-tmin+1
etr = t2-tmin+1
ita = t3-tmin+1
eta = t4-tmin+1


date_i, date_f, time = searchtime_pos(tmin-1, tmax-1)
##print('date_i={}, date_f={}, time={}'.format(date_i, date_f, time))
date1 = date_i
date2 = date1[:5]+iadate
# ##print('date2 => ', date2)
# ##print('time => ', time)
varctl_idx=''


path= _workdir+'/outputs/prueba/'
if not os.path.exists(path + ip ):
    os.makedirs(path+ip)
    f = open(path+ip+'/file_query.log', 'w')
    f.close()

path = path+ip


def save_log():
    with open(path+'/file_query.log', 'r') as old_query:
        data = old_query.readlines()
        old_query.close()
    with open(path+'/file_query.log', 'w') as new_query:
        for line in data:
            new_query.write(line)
        # new_query.write('\n')
        new_query.write(nowdate+'\n')
        new_query.write(server_cad+'\n')
        new_query.write('\n')
        new_query.close()


save_log()

# ofile1 = 'mando'

option, average = operation(tmin, tmax, '/'.join([path, ofile1]))
# #print('option: {}\naverage: {}'.format(option, average))
# exit()


write_ctl(path, ofile1,nx,xfirst,xinc,ny,yfirst,yinc,time,date1,undef,var1,varctl_idx)
##print()
# ifile1='/'.join([path,ofile1+'.grb'])
# if(dindex=='spi'):
#     ifile2='/'.join([path,ofile1+'.grb'])
# else:
#     ifile2='/'.join([path,ofile1+'.grb'])

ifile1='/'.join([path,ofile1+'.grb'])
if(dindex=='spi'):
    ifile2='/'.join([path,ofile1+'.grb'])
else:
    ifile2='/'.join([path,ofile1+'.grb'])




##print('ifile1: {}'.format(ifile1))
##print('ifile2: {}'.format(ifile2))
##print()
totlen=time
tlen=etr-itr+1
tlen2=eta-ita+1

##print('tlen:{}, tlen2:{}'.format(tlen, tlen2))

ofile='_'.join([modeldir,str(ilon),str(ilat),str(elon),str(elat),str(iryear)+'-'+str(eryear),str(iayear)+'-'+str(eayear),option,dindex+str(spilen),'usrthreshold',str(spilim),'nmonths',str(evenlim),'thresholdcategories'+str(extdrou),str(sevdrou),str(moddrou),'rawdataindex-output'])
pfile='_'.join([modeldir,str(ilon),str(ilat),str(elon),str(elat),str(iryear)+'-'+str(eryear),str(iayear)+'-'+str(eayear),option,dindex+str(spilen),'usrthreshold',str(spilim),'nmonths',str(evenlim),'thresholdcategories'+str(extdrou),str(sevdrou),str(moddrou),'droughtstats-output'])


ovar=dindex

# #print('ofile: {}'.format(ofile))
# #print('pfile: {}'.format(pfile))
# #print('ovar: ',ovar)
# #print()

def return_args(lista):

    lista_args = []
    for i,v in enumerate(lista):
        if not isinstance(v,str):
            lista_args.append(str(v))
        else:
            lista_args.append(v)

    return lista_args



cmdfiles=return_args([ifile1,ifile2,'/'.join([path,ofile+'.grb']),'/'.join([path,pfile])])
cmdlinepar= return_args([totlen,tlen,itr,etr,tlen2,ita,eta,'12',nx,ny,yfirst,yinc,undef])
cmdlinelim = return_args([spilen,evenlim,spilim])
cmdlinedrou=return_args([extdrou,sevdrou,moddrou])
cmdline = cmdfiles + cmdlinepar + cmdlinelim + cmdlinedrou
# #print(cmdline)
# exit()

if (dindex=='spi'):
    # Compile FORTRAN
    # #print('aqui esta el spi.f90', '/'.join([_workdir, progdir,'spi.f90']))
    # exit()
    subprocess.run('/usr/bin/gfortran {} -o {}'.format('/'.join([_workdir, progdir,'spi.f90']),'/'.join([_workdir, progdir,'spi'])), shell=True)
    os.chdir(home)
    # # # Run FORTRAN (spi)
    spif90 = subprocess.run('./data4drought/programs/spi {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'
                .format(
                        cmdline[0],
                        cmdline[1],
                        cmdline[2],
                        cmdline[3],
                        cmdline[4],
                        cmdline[5],
                        cmdline[6],
                        cmdline[7],
                        cmdline[8],
                        cmdline[9],
                        cmdline[10],
                        cmdline[11],
                        cmdline[12],
                        cmdline[13],
                        cmdline[14],
                        cmdline[15],
                        cmdline[16],
                        cmdline[17],
                        cmdline[18],
                        cmdline[19],
                        cmdline[20],
                        cmdline[21],
                        cmdline[22],

                            ), shell=True, capture_output=True)

    # #print('spif90.stderr:\n',spif90.stderr)
    # #print()
    # #print('spif90.stdout:\n',spif90.stdout)

    # #print('runned spi FORTRAN script')


else:
    # Compile FORTRAN
    subprocess.run('/usr/bin/gfortran  {} -o {}'.format('/'.join([_workdir, progdir,'spei2.f90']), '/'.join([_workdir, progdir,'spei2'])), shell=True)
    os.chdir(home)
    # # # Run FORTRAN (spei)
    speif90 = subprocess.run('./data4drought/programs/spei2 {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}'
                .format(
                        cmdline[0],
                        cmdline[1],
                        cmdline[2],
                        cmdline[3],
                        cmdline[4],
                        cmdline[5],
                        cmdline[6],
                        cmdline[7],
                        cmdline[8],
                        cmdline[9],
                        cmdline[10],
                        cmdline[11],
                        cmdline[12],
                        cmdline[13],
                        cmdline[14],
                        cmdline[15],
                        cmdline[16],
                        cmdline[17],
                        cmdline[18],
                        cmdline[19],
                        cmdline[20],
                        cmdline[21],
                        cmdline[22],

                            ), shell=True, capture_output=True)

    # #print(speif90.stdout)
    # #print('runned spei FORTRAN script')


# ##print()
# remove = subprocess.run('rm {}'.format('/'.join([path,'*.nc'])), shell=True, capture_output=True) 
# # remove = subprocess.run('rm {}'.format('/'.join([_workdir,'outputs','prueba','*.nc'])), shell=True, capture_output=True) 
# # #print('old nc eliminados')
# ##print()




write_ctl(path,ofile,nx,xfirst,xinc,ny,yfirst,yinc,tlen2,date2,undef,ovar,spilen)
# #print('write_ctl {}'.format(ofile))
# exit()

pfile1=pfile+'_1'
ovar='totmon'
rc=write_ctl(path,pfile1,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar,spilen)

pfile1=pfile+'_2'
ovar='maxmon'
rc=write_ctl(path,pfile1,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar,spilen)

pfile1=pfile+'_3'
ovar='nevent'
rc=write_ctl(path,pfile1,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar,spilen)

pfile1=pfile+'_4'
ovar='prcmon'
rc=write_ctl(path,pfile1,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar,spilen)




pfile2=pfile+'d_1'
ovar2='extdry'
rc=write_ctl(path,pfile2,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar2,spilen)

pfile2=pfile+'d_2'
ovar2='sevdry'
rc=write_ctl(path,pfile2,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar2,spilen)

pfile2=pfile+'d_3'
ovar2='moddry'
rc=write_ctl(path,pfile2,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar2,spilen)

pfile2=pfile+'d_4'
ovar2='normal'
rc=write_ctl(path,pfile2,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar2,spilen)

pfile2=pfile+'d_5'
ovar2='modwet'
rc=write_ctl(path,pfile2,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar2,spilen)

pfile2=pfile+'d_6'
ovar2='sevwet'
rc=write_ctl(path,pfile2,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar2,spilen)

pfile2=pfile+'d_7'
ovar2='extwet'
rc=write_ctl(path,pfile2,nx,xfirst,xinc,ny,yfirst,yinc,1,date2,undef,ovar2,spilen)



##print('yaaaaaaaaaaaaaaaaaaaa CTL')


ofile_ctl='/'.join([path,ofile])
ofile_nc='/'.join([path,ofile])
pfile = '/'.join([path,pfile])



##print('ofile_ctl: {}'.format(ofile_ctl))
##print('ofile_nc: {}'.format(ofile_nc))
##print('pfile: {}'.format(pfile))


subprocess.run('cdo -s -f nc import_binary {}.ctl {}.nc'.format(ofile_ctl, ofile_nc), shell=True)

# #print('{}.ctl to {}.nc'.format(ofile_ctl, ofile_nc))
# exit()

for i in range(1,5):
    subprocess.run('cdo -s -f nc import_binary {}_{}.ctl {}_{}.nc'.format(pfile,i, pfile,i), shell=True)


subprocess.run('cdo -s merge {}_1.nc {}_4.nc {}_2.nc {}_3.nc {}.nc'.format(pfile, pfile, pfile, pfile, pfile), shell=True)


for i in range(1,8):
    subprocess.run('cdo -s -f nc import_binary {}d_{}.ctl {}d_{}.nc'.format(pfile,i, pfile,i), shell=True)


subprocess.run('cdo -s merge {}d_1.nc {}d_2.nc {}d_3.nc {}d_4.nc {}d_5.nc {}d_6.nc {}d_7.nc {}kind_drought.nc'.format(pfile, pfile, pfile, pfile, pfile, pfile, pfile, pfile), shell=True)


#  se borran todos los files grb, ctl y los 4 netcdf temporales
subprocess.run('rm {}'.format(path+'/*.grb'), shell=True)
subprocess.run('rm {}'.format(path+'/*.ctl'), shell=True)
for i in range(1,5):
    subprocess.run('rm {}_{}.nc'.format(pfile,i), shell=True)
for i in range(1,8):
    subprocess.run('rm {}d_{}.nc'.format(pfile,i), shell=True)



# # definiendo nombre de los 3 files de salida
# ofile_raw = ofile_nc
# ofile_count = pfile
# ofile_kind = pfile+'kind_drought'





# subprocess.run('/usr/bin/grads -lbc "run scripts.gs {} {} {} {} {} {} {} {} {} {} {}" >> log.txt'
#                 .format('adrian', 30, avg, ilat, elat, ilon, elon, irdate, erdate, iadate, eadate), 
#                 shell=True, encoding='utf-8')
# ##print('runned grads script')



# os.environ['GADDIR'] = '/usr/lib/cgi-bin/grads-2.1.a2.oga.1/Classic/data'
# os.environ['GASCRP'] = '/usr/lib/cgi-bin/grads-2.1.a2.oga.1/Classic/scripts'
# os.environ['GAUDFT'] = '/usr/lib/cgi-bin/grads-2.1.a2.oga.1/Classic/data'

# subprocess.run('/usr/lib/cgi-bin/grads-2.1.a2.oga.1/Classic/bin/grads -lbc "run scripts.gs {} {} {} {} {} {} {} {} {} {} {}" >> log.txt'
#                 .format('adrian', 30, avg, ilat, elat, ilon, elon, irdate, erdate, iadate, eadate), 
#                 shell=True, encoding='utf-8')



# minlat_ref, maxlat_ref = 5, 29
# minlon_ref, maxlon_ref = -94.5, -50

# minlat_ref, maxlat_ref = ilat, elat 
# minlon_ref, maxlon_ref = ilon, elon

# pons_ini, pos_end = 13, 493
# pons_ini, pos_end = 40, 90

# precip_py = save_mean(var=var1, ini=tmin, end=tmax, minlat=ilat, minlon=ilon, maxlat=elat, maxlon=elon)

# precip_grb = read_from_grb('./data4drought/outputs/prueba/ficherodeprueba.grb')

# file_precip_txt = open('./data4drought/outputs/prueba/output_fortran.txt', 'r')
# raw_precip_txt = file_precip_txt.read().splitlines()
# precip_txt = []
# for value in raw_precip_txt:
#     precip_txt.append(float(value))


# ##print(precip_py.shape)
# dif_grb_py = precip_grb-precip_py

# for (grb,py,diff) in zip(precip_grb, precip_py, dif_grb_py):
#     ##print('grb:{}, py:{}, grb-py:{}'.format(grb,py,diff))

# ##print()
# ##print(dif_grb_py.min(), dif_grb_py.max())
    
# import matplotlib.pyplot as plt

# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))

# ax1.plot(precip_grb,'-',label='grb')
# ax1.plot(precip_py,'-', label='py')
# ax1.legend()
# ax2.plot(dif_grb_py, '-', label='bias')
# ax2.legend()
# plt.ylim([-10,10])
# plt.savefig('./data4drought/outputs/prueba/averages.png',dpi=300)

# ##print('grb: {} py: {} txt: {}'.format(len(precip_grb), len(precip_py), len(precip_txt)))





# out=subprocess.run('tar -czvf uploads/AAA.tar.gz uploads/*.png', shell=True, capture_output=True)

# dir_nc_files = '/'.join([_workdir,'outputs', 'prueba'])
os.chdir(path)
# current_dir = os.getcwd()
# #print(dir_nc_files)
# #print(current_dir)
files = glob.glob('*.nc')



nowdate += '.tar.gz'
namefile = ip+'-'+'ncfiles'+'_'+nowdate

subprocess.run('tar -czf {} *.nc'.format(namefile), shell=True)


subprocess.run('mv {} {}'.format(namefile,root_dir+'/uploads/targz_files'), shell=True)

subprocess.run('rm *.nc', shell=True)

print('adrian')
# print(namefile)
