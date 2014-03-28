for archivo in *.data
do
    echo "mirando $archivo"
    salida=$(awk -F" " '{print $2}' $archivo | grep "1\.00" | head -n 0)
    if test -z $salida
	then
	echo "es ejecutable"
    fi
done