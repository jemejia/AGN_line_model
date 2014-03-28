



for nombre in $(cut -f 2 ../quasar_data.txt)
do
    
    mkdir $nombre
    seleccion=$(ls *.png | grep $nombre)  
    mv $seleccion $nombre
    #salida=$(awk -F" " '{print $2}' $archivo | grep "1\.00" | head -n 0)
    #if test -z $salida
    #then
    #echo "es ejecutable"
    #fi
done