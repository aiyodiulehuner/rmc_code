for i in 80 60
do
echo $i
../cofirank/dist/cofirank-deploy cofirank_ns_$i.config
done
