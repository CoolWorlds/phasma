#!/bin/bash


targetlist='./targetlist.txt'
istart=2
iend=`wc -l $targetlist`
iend=`echo $iend | cut -f1 -d ' '`
zero='0'
rawLC='rawLC'

#download Kepler fits files for targets in targetlist.txt
for i in `seq $istart $iend`; do
  # ===================================
  # ========= FILE DEFINITION =========
  # ===================================
  KOIname=`head -$i $targetlist | tail -1 | cut -f1 -d$'\t'` #tab-delimited
  KIC=`head -$i $targetlist | tail -1 | cut -f2 -d$'\t'` #tab-delimited
  
  echo ' '
  echo 'PROCESSING '$KOIname' ('$(($i - 1))' of '$(($iend - 1))')'
  # Define a KIC name which has been zero padded
  KICchars=`echo $KIC | wc -c`
  KICchars=$(($KICchars - 1))
  zeros=$((9 - $KICchars))
  if [ "$zeros" -eq 0 ]; then
    KICname=$KIC
  fi
  if [ "$zeros" -eq 1 ]; then
    KICname=`echo $zero$KIC`
  fi
  if [ "$zeros" -eq 2 ]; then
    KICname=`echo $zero$zero$KIC`
  fi
  if [ "$zeros" -eq 3 ]; then
    KICname=`echo $zero$zero$zero$KIC`
  fi
  if [ "$zeros" -eq 4 ]; then
    KICname=`echo $zero$zero$zero$zero$KIC`
  fi
  if [ "$zeros" -eq 5 ]; then
    KICname=`echo $zero$zero$zero$zero$zero$KIC`
  fi
  if [ "$zeros" -eq 6 ]; then
    KICname=`echo $zero$zero$zero$zero$zero$zero$KIC`
  fi
  if [ "$zeros" -eq 7 ]; then
    KICname=`echo $zero$zero$zero$zero$zero$zero$zero$KIC`
  fi
  if [ "$zeros" -eq 8 ]; then
    KICname=`echo $zero$zero$zero$zero$zero$zero$zero$zero$KIC`
  fi
  KICshort=${KICname:0:4}
  # Check if a directory has already been made or not
  folder=`echo "./KOIs/"$KOIname`
  dirfound=`ls | grep -ci $folder`
  if [ "$dirfound" -eq 0 ]; then
    echo 'Making a new directory '$folder
    mkdir $folder
  else
    echo 'No directory creation needed'
  fi
  # Check if the directory has a rawLC subdirectory
  cd $folder
  dirfound=`ls | grep -ci $rawLC`
  if [ "$dirfound" -eq 0 ]; then
    mkdir $rawLC
  fi
  cd ../..
  pwd
  # ===================================
  # ========== FILE DOWNLOAD ==========
  # ===================================
  # Check if files have been downloaded yet or not
  nfits=`ls $folder"/"$rawLC"/" | grep -c fits`
  if [ "$nfits" -eq 0 ]; then
    cd $folder
    cd $rawLC
    ftpfolder='http://archive.stsci.edu/pub/kepler/lightcurves//'
    ftpfolder=`echo $ftpfolder$KICshort'/'$KICname'/'`
    echo 'Downloading files for '$folder
    wget -q -nH --cut-dirs=6 -r -l0 -c -N -np -R 'index*' -erobots=off $ftpfolder
    cd ../../..
    pwd
    nfits=`ls $folder"/"$rawLC"/" | grep -c fits`
    echo 'Downloaded '$nfits' new files'
  else
    echo 'No need to download. We already have '$nfits' files'
  fi
done
echo 'Done'

#get data out of fits files and run phasma detrending
for i in `seq $istart $iend`; do
  KOIname=`head -$i $targetlist | tail -1 | cut -f1 -d$'\t'` #tab-delimited
  KIC=`head -$i $targetlist | tail -1 | cut -f2 -d$'\t'` #tab-delimited

  echo ' '
  echo 'PROCESSING '$KOIname' ('$(($i - 1))' of '$(($iend - 1))')'
  
  # Check if the planet has already been processed (i.e. if "intermediate" and "final" files already exist)
  int_filename=`echo "./KOIs/"$KOIname"/"$KOIname"_phasma_intermediate.txt"`
  final_filename=`echo "./KOIs/"$KOIname"/"$KOIname"_phasma_final.csv"`
  
  if [ -e "$int_filename" ]; then
    echo 'No need to generate intermediate file for '$KOIname'. File exists'

  else
    echo 'Generating intermediate file for '$KOIname
    
    while read line; do
      lineKIC=`echo $line | cut -f2 -d ','`
      lineKOI=`echo $line | cut -f3 -d ','`
      
      if ([ "$lineKIC" == "$KIC" ] && [ "$lineKOI" == "$KOIname" ]); then
        NEAdisp=`echo $line | cut -f5 -d ','`
        echo $NEAdisp
        Kepdisp=`echo $line | cut -f8 -d ','`
        echo $Kepdisp
        t_orbit=`echo $line | cut -f16 -d ','`
        echo $t_orbit
        t_midpt=`echo $line | cut -f17 -d ','`
        echo $t_midpt
        t_14=`echo $line | cut -f22 -d ','`
        echo $t_14
        rp_rstar=`echo $line | cut -f25 -d ','`
        echo $rp_rstar
        inc=`echo $line | cut -f30 -d ','`
        echo $inc
        nplanet=`echo $line | cut -f43 -d ','`
        echo $nplanet
        DVlink=`echo $line | cut -f53 -d ','`
        echo $DVlink
        stellar_teff=`echo $line | cut -f54 -d ','`
        echo $stellar_teff
        stellar_rad=`echo $line | cut -f57 -d ','`
        echo $stellar_rad
        stellar_mass=`echo $line | cut -f58 -d ','`
        echo $stellar_mass
        break
      fi
    done <KOIarchive_all_heartbeat.csv

    folder=`echo "./KOIs/"$KOIname`
    LCfolder=`echo "./KOIs/"$KOIname"/rawLC"`

    intermediatefile=`echo "./KOIs/"$KOIname"/"$KOIname"_phasma_intermediate.txt"`
    finalfile=`echo "./KOIs/"$KOIname"/"$KOIname"_phasma_final.csv"`
    binnedfile=`echo "./"$KOIname".csv"`

    for f in $LCfolder/*.fits; do
      #python extract_fits.py $f $intermediatefile $t_orbit $t_14 $t_midpt 
      python phasma.py $f $intermediatefile $t_orbit $t_14 $t_midpt 
    done
  fi

  if [ -e "$final_filename" ]; then
    echo 'No need to generate final file for '$KOIname'. File exists'

  else
    echo 'Generating final file for '$KOIname
    
    while read line; do
      lineKIC=`echo $line | cut -f2 -d ','`
      lineKOI=`echo $line | cut -f3 -d ','`
      
      if ([ "$lineKIC" == "$KIC" ] && [ "$lineKOI" == "$KOIname" ]); then
        NEAdisp=`echo $line | cut -f5 -d ','`
        echo $NEAdisp
        Kepdisp=`echo $line | cut -f8 -d ','`
        echo $Kepdisp
        t_orbit=`echo $line | cut -f16 -d ','`
        echo $t_orbit
        t_midpt=`echo $line | cut -f17 -d ','`
        echo $t_midpt
        t_14=`echo $line | cut -f22 -d ','`
        echo $t_14
        rp_rstar=`echo $line | cut -f25 -d ','`
        echo $rp_rstar
        inc=`echo $line | cut -f30 -d ','`
        echo $inc
        nplanet=`echo $line | cut -f43 -d ','`
        echo $nplanet
        DVlink=`echo $line | cut -f53 -d ','`
        echo $DVlink
        stellar_teff=`echo $line | cut -f54 -d ','`
        echo $stellar_teff
        stellar_rad=`echo $line | cut -f57 -d ','`
        echo $stellar_rad
        stellar_mass=`echo $line | cut -f58 -d ','`
        echo $stellar_mass
        break
      fi
    done <KOIarchive_all_heartbeat.csv

    folder=`echo "./KOIs/"$KOIname`
    LCfolder=`echo "./KOIs/"$KOIname"/rawLC"`

    intermediatefile=`echo "./KOIs/"$KOIname"/"$KOIname"_phasma_intermediate.txt"`
    finalfile=`echo "./KOIs/"$KOIname"/"$KOIname"_phasma_final.csv"`
    binnedfile=`echo "./"$KOIname".csv"`

    python phasma_masterfile_processing.py $intermediatefile $finalfile $binnedfile $t_orbit $t_14 $t_midpt $stellar_rad $stellar_mass $stellar_teff $rp_rstar $inc $nplanet $NEAdisp $Kepdisp $DVlink
  fi
done

echo 'Done done done'

