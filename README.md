C implementation of merge sort, radix sort and parallel radix sort on 32 bit unsigned integers

![](https://cdn.rawgit.com/983/RadixSort/master/sorted_data.svg "")
![](https://cdn.rawgit.com/983/RadixSort/master/reversed_data.svg "")

Performance of qsort and merge sort degrades as data becomes more unordered,
while running time of radix sort is much less affected.

![](https://cdn.rawgit.com/983/RadixSort/master/90_percent_sorted_data.svg "")
![](https://cdn.rawgit.com/983/RadixSort/master/random_data.svg "")

Parallel radix sort is only worth it if the number of elements is large.

![](https://cdn.rawgit.com/983/RadixSort/master/random_data_2.svg "")

Tested on an i5-4210U with 2 cores. Algorithm will scale gracefully with higher number of cores.

For very small arrays, qsort, merge sort and sorting networks should be favoured over radix sort.
