rm tmp.txt
if [ ! -d viruses_bins_functional_marker/ ];then mkdir viruses_bins_functional_marker;fi
cat $1 |sed -n "/${2}/p"|cut -f 2|while read line;do 
j=${line%__*};
echo -e $j'\t'$line >> tmp.txt;
done
cat tmp.txt |uniq > bins_functional_marker/${2}_${3}.txt

