TARGET = bf.out linked_list.out linked_list_simd1x4.out linked_list_simd4x1.out

all: $(TARGET)

bf.out: force.cpp
	icpc -march=core-avx2 -std=c++11 -DBRUTE_FORCE -O3 $< -o $@

linked_list.out: force.cpp
	icpc -march=core-avx2 -std=c++11 -O3 $< -o $@

linked_list_simd1x4.out: force.cpp
	icpc -march=core-avx2 -std=c++11 -DSIMD1x4 -O3 $< -o $@

linked_list_simd4x1.out: force.cpp
	icpc -march=core-avx2 -std=c++11 -DSIMD4x1 -O3 $< -o $@

clean:
	rm -f $(TARGET)

test: $(TARGET)
	./linked_list.out 2> linked_list.dat
	./linked_list_simd1x4.out 2> linked_list_simd1x4.dat
	./linked_list_simd4x1.out 2> linked_list_simd4x1.dat
	diff linked_list.dat linked_list_simd1x4.dat
	diff linked_list_simd1x4.dat linked_list_simd4x1.dat
