# %%
# f'GPH anomaly from the {syear}-{eyear} at {prs}hPa normalized by the standard deviation'

import calendar
from matplotlib import axis, cm
import numpy as np
import math
from datetime import date
import os
import matplotlib.pyplot as plt

# ==========初期値============
startMonth = 9
endMonth = 11
syear = 2010
eyear = 2020
name = 'hgt'
prsInd = 31
meanLat = [-60,-50]

# ===========定数================
pcord = np.array([1000,975,950,925,900,875,850,825,800,775,750,700,
        650,600,550,500,450,400,350,300,250,225,200,175,150,125,100,70,
        50,30,20,10,7,5,3,2,1])
kind = 'std'
grid = 1.25
phicord = np.arange(-90,91,grid)*(math.pi/180.)
xcord = np.arange(-180,181,grid)
prs = pcord[prsInd]


def character(meanlatrange):
    characterMLR = []
    for lat in meanlatrange:
        if lat > 0:
            characterMLR.append(f'{abs(lat)}N')
        else:
            characterMLR.append(f'{abs(lat)}S')
    return characterMLR

def JRA_latInd(meanLat):
    meanLat = np.array(meanLat)
    meanLatInd = np.array(meanLat/grid + (90/grid), np.int32)
    yarray = np.arange(meanLat[0],meanLat[1]+grid,grid)
    phiarray = np.radians(yarray)
    cosphi = np.cos(phiarray)
    return meanLatInd,cosphi


def check_uruYear(year):
    if year % 4 ==0:
        allcday = 366
    else:
        allcday = 365
    return allcday

def dayIndRange(year,sMonth,eMonth):
    sdayNumInd = (date(year,sMonth,1)-date(year,1,1)).days
    edayNumInd = (date(year,eMonth+1,1)-date(year,1,1)).days-1
    return sdayNumInd, edayNumInd

def exAnyPrsAndDate_data(year,name,prsInd,sdayNumInd,edayNumInd):
    allcday = check_uruYear(year)
    meanLatInd,cosphi = JRA_latInd(meanLat)
    file = f'D:/data/JRA55/{name}/anl_p_{name}.{year}.bin'
    f = open(file, 'rb')
    print(f'{name}/{year} 読み込み中')
    allArray = np.fromfile(f,dtype='>f').reshape(allcday,37,145,288)
    print(f'{name}/{year}読み込み官僚')
    anyArray = allArray[sdayNumInd:edayNumInd+1,prsInd,meanLatInd[0]:meanLatInd[1]+1,:]
    anyArray = np.average(anyArray, axis=1, weights=cosphi) #anyArray.shape=(anydaynum,288)
    daycount = 0
    for mon in range(startMonth,endMonth+1):
        dayNum = calendar.monthrange(year,mon)[1]
        lonArray = np.mean(anyArray[daycount:daycount+dayNum],axis=0)
        if mon == startMonth:
            monthArray_2d = lonArray[np.newaxis]
        else:
            monthArray_2d = np.append(monthArray_2d,lonArray[np.newaxis],axis=0)
        daycount += dayNum

    return monthArray_2d

def checkDir(filename):
    dirList = filename.split('/')
    for dir in dirList:
        if dirList.index(dir) == 0:
            pathname = dir
        elif not dir == dirList[-1]:
            pathname += f'/{dir}'
            if not os.path.exists(pathname):
                os.mkdir(pathname)
                print(f'make directory "{pathname}"')

def calStDev(name,syear,eyear,savename):
    for year in range(syear,eyear+1):
        sdayNumInd, edayNumInd = dayIndRange(year,startMonth,endMonth)
        array_2d = exAnyPrsAndDate_data(year,name,prsInd,sdayNumInd,edayNumInd)
        if year == syear:
            array_3d = array_2d[np.newaxis]
        else:
            array_3d = np.append(array_3d,array_2d[np.newaxis],axis=0)
    mean_2d = np.mean(array_3d,axis=0)   
    std_2d = np.std(array_3d,axis=0,ddof=1)
    normalized_by_Std = (array_3d-mean_2d)/std_2d
    
    for i in range(eyear-syear+1):
        if i == 0:
            normalized_by_Std_2d = normalized_by_Std[i]
        else:
            normalized_by_Std_2d = np.append(normalized_by_Std_2d,normalized_by_Std[i],axis=0)
    checkDir(savename)
    np.save(savename,normalized_by_Std_2d)
    return normalized_by_Std_2d

def changeJRA55Lon(array_2d): #経度座標を0to358.75から-180to180に変える
    half_east = array_2d[:,:73]
    half_west = array_2d[:,72:]
    newarray_2d = np.append(half_west,half_east,axis=1)
    return newarray_2d

def yNum():
    allMonthNum = (eyear-syear+1)*(endMonth-startMonth+1)
    monthNum = (endMonth-startMonth+1)
    return allMonthNum, monthNum

def convertValue(value):
    if value < -4:
        ind = 0
    elif value >= 4:
        ind = -1
    else:
        ind = int((value+4)*2)
    return ind


def draw(array_2d,title,allMonthNum, monthNum,picturename,alllon):
    # カラーバーを作るための適当な図
    fig2,ax2 = plt.subplots(figsize=(1,1))
    cm16 = plt.get_cmap('jet',16)
    array = np.arange(20).reshape(4,5)
    im = plt.imshow(array,cmap=cm16,vmin=-4,vmax=4)

    # 本番図の描写開始
    fig, ax = plt.subplots(figsize=(6,6),facecolor='#fff')

    # y軸---------------
    yday = np.arange(0,allMonthNum+1,monthNum)
    ydayStr = np.array(np.arange(syear,eyear+2),dtype=np.str_)
    ax.set_ylim(0,allMonthNum)
    ax.set_yticks(yday)
    ax.set_yticklabels(ydayStr)
    ax.set_ylabel(f'year ( {startMonth} - {endMonth} )')

    # x軸---------------
    xlon = np.arange(0,290,48)
    xlonStr =  np.array(np.arange(-180,181,60),dtype=np.str_)
    ax.set_xlim(0,289)
    ax.set_xticks(xlon)
    ax.set_xticklabels(xlonStr)
    ax.set_xlabel('lon')

    # カラーバー設置-----------
    # axpos = ax.get_position()
    # cbar_ax = fig.add_axes([0.81, axpos.y0, 0.02, axpos.height])
    # plt.colorbar(im,cax=cbar_ax)
    plt.colorbar(im)
    
    # 図の配置----------
    plt.subplots_adjust(right=0.78)
    
    # タイトル----------
    plt.title(title)

    # データ描画----------------
    for month in range(allMonthNum):
        for lon in range(alllon):
            anycolor = cm16(convertValue(array_2d[month,lon]))
            if lon == 0:
                ax.axvspan(
                    0,
                    0.5,
                    0.2,
                    0.6,
                    color = anycolor
                    )
            elif lon == alllon:
                ax.axvspan(
                    lon-0.5,
                    lon,
                    month / allMonthNum,
                    (month + 1) / allMonthNum,
                    color = anycolor
                    )
            else:
                ax.axvspan(
                    lon - 0.5,
                    lon + 0.5,
                    month / allMonthNum,
                    (month + 1) / allMonthNum,
                    color = anycolor,
                    )

    # x軸平衡の実線と破線-----------
    for i in range(1,34):
        if i % 3 == 0:
            lineKind = "solid"
            linewidth = 0.9
            color = "black"
        else:
            lineKind = "dashed"
            linewidth = 0.5
            color = 'grey'
        plt.hlines(i,0,289,
            color=color, 
            linewidth=linewidth, 
            linestyles=lineKind,
            )
    
    # 図の保存------------------           
    # if not os.path.exists(picturename):
    checkDir(picturename)
    plt.savefig(picturename)
    
    # 図のディスプレイ表示-------------
    plt.show()
    

def main():
    allMonthNum, monthNum = yNum()
    characterMLR = character(meanLat)
    savename = f'D:/data/JRA55/{name}/{syear}-{eyear}/{syear}-{eyear}_{characterMLR[0]}-{characterMLR[1]}_{prs}hPa_{name}_{kind}.npy'
    
    if not os.path.exists(savename):
        std_2d =  calStDev(name,syear,eyear,savename)
    else:
        std_2d = np.load(savename)
    new_nbsdt_2d = changeJRA55Lon(std_2d)

    picturetitle = f'Normalized GPH {characterMLR[0]}-{characterMLR[1]} {prs}hPa'
    picturename = f'D:/picture/study/JRA55/blocking/Normalized_GPH_{characterMLR[0]}-{characterMLR[1]}_{prs}hPa'
    draw(new_nbsdt_2d,picturetitle,allMonthNum,monthNum,picturename,288)

if __name__ == '__main__':
    main()
     