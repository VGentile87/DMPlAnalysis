#!bin/bash
echo -e "Type: \n 1 File with images \n 2 File without images"
read var_im
echo Hai scelto $var_im
FILE_IM='dm_tracks_cl.dm.root'
FILE_noIM='dm_tracks.dm.root'

if [ $var_im = '1' ] ; then
    echo "ciao"
    file_name=$FILE_IM

elif [ $var_im = '2' ] ; then
    file_name=$FILE_noIM
fi

echo Il tuo file di analisi è $file_name

if [ ! -f $file_name ]; then
    echo "Il file non esiste"

else

    echo -e "Type:\n 1 for MakeClass + Analysis \n 2 for Only Analysis \n 3 for GIF macro"
    read varname
    echo Hai scelto $varname
    if [ $varname = '1' ]; then
	cp  $DMPLA/start_v2.C .
	root -l $file_name <<EOC
.x start_v2.C()
.q
EOC
	rm AutoDict_vector_vector_string_allocator_string_____*
	sed -n "36, 40p" myData.h > tmp_file
	sed -i 27rtmp_file myData_v6.C
	if [ $var_im = '2' ]; then
	    #echo "QUI CI ARRIVO"
	    sed -i '0,/images/! {0,/images/ s=images=//images=}' myData_v6.C
	fi
	root -l $file_name <<EOC
gSystem->Load("libDMRoot");
gROOT->ProcessLine(".L myData_v6.C++");
myrun();
.q
EOC
	echo "Look histograms in data.root file and check if settings.mac is ok"
	more settings.mac
	echo "If some setting values are wrong, please change them in settings.mac of this directory"
	echo "Then type 'Only Analysis' option to run again the analysis tool"
    elif [ $varname =  '2' ]; then
	cp $DMPLA/macros/myData_v6.C .
	sed -n "36, 40p" myData.h > tmp_file
	sed -i 27rtmp_file myData_v6.C
	if [ $var_im = '2' ]; then
	    echo "QUI CI ARRIVO"
	    sed -i '0,/images/! {0,/images/ s=images=//images=}' myData_v6.C
	fi
	#sed -i '0,/images/! {0,/images/ s=images=//images=}' myData_v6.C
	root -l $file_name <<EOC
gSystem->Load("libDMRoot");
gROOT->ProcessLine(".L myData_v6.C++");
myrun();
.q
EOC
	more settings.mac
    elif [ $varname =  '3' ]; then
	cp $DMPLA/macros/run_gif.C .
	cp $DMPLA/macros/dgr.C .
	echo 'dgr.C and run_gif.C copied in this directory'
    
    else
	echo "Number incorrect...bye"
    fi
fi
