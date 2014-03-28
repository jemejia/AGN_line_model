for empty in $(du -h | grep '0B' | cut -f 2)
do
    rm -rf $empty
done