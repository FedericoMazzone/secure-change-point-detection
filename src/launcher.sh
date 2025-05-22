# for epsilon in 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0
# repeating each experiment 10 times
# run ./build/clustering -f data/generated/points_d2_n10000_k8.csv -k 8 -r 2 -e epsilon -t
# store the output and error of each run in a file

for epsilon in 0.5 1.0 2.0 3.0 4.0 5.0
do
    for i in $(seq 1 30)
    do
        echo "Running epsilon=$epsilon, iteration=$i..."
        filename="s1_epsilon${epsilon}_run${i}.log"
        ./build/clustering -f data/s1.csv -k 15 -r 2 -e $epsilon -t -id 5 > $filename 2>&1
    done
done
