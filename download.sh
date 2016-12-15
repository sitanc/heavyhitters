#!/bin/bash
# time=130000
# curl_address="https://data.caida.org/datasets/passive-2016/equinix-chicago/20160121-130000.UTC/equinix-chicago.dirA.20160121-125911.UTC.anon.pcap.gz"
# dir_path="data/$time"
# pcapgz_path="$dir_path/data.pcap.gz"
# pcap_path="$dir_path/data.pcap"
# output_path="$dir_path/output.txt"
# mkdir $dir_path
# curl $curl_address -o $pcapgz_path -u "sitanc@mit.edu:Lotrfreak2131995"
# gunzip -f $pcapgz_path
# tcpdump -l -r $pcap_path -n ip | awk '{ print gensub(/(.*)\..*/,"\\1","g",$3), $4, gensub(/(.*)\..*/,"\\1","g",$5) }' > $output_path
# rm $pcap_path

for time in {130000,130100,130400,131600,133200,134200,135200,140200}
do
    curl_address="https://data.caida.org/datasets/passive-2016/equinix-chicago/20160121-130000.UTC/equinix-chicago.dirA.20160121-$time.UTC.anon.pcap.gz"
    dir_path="data/$time"
    pcapgz_path="$dir_path/data.pcap.gz"
    pcap_path="$dir_path/data.pcap"
    output_path="$dir_path/output_longer.txt"
    mkdir $dir_path
    curl $curl_address -o $pcapgz_path -u "sitanc@mit.edu:Lotrfreak2131995"
    gunzip -f $pcapgz_path
	tcpdump -l -r $pcap_path -n ip | awk '{ print gensub(/(.*)\..*/,"\\1","g",$3), $4, gensub(/(.*)\..*/,"\\1","g",$5) }' > $output_path
    rm $pcap_path
done

python hist.py > final_data.txt