for i in 20 40 60 80
do
echo $i
../cofirank/dist/cofirank-deploy cofirank_ns_$i.config
done
