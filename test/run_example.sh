#! \bin\bash

# Make fure it has been compiled before using it
fpath="../src"

for i in {1..10}

do
# FOR seq version
${fpath}/TEST './Table/dna3_40_4_16_new.table' './Input/dna/small.in' 40 4 8 8 16 2 0 >> RECORD

# FOR Multi-thread Seq 
${fpath}/TEST './Table/dna3_40_4_16_new.table' './Input/dna/small.in' 40 4 8 8 16 2 1 >> RECORD

# FOR spec nonSIMD version
${fpath}/TEST './Table/dna3_40_4_16_new.table' './Input/dna/small.in' 40 4 8 8 16 2 2 >> RECORD

# FOR spec AVX version
${fpath}/TEST './Table/dna3_40_4_16_new.table' './Input/dna/small.in' 40 4 8 8 16 2 3 >> RECORD

# FOR sse with range coalescing version
${fpath}/TEST './Table/dna3_40_4_16_new.table' './Input/dna/small.in' 40 4 8 8 16 2 6 >> RECORD

# FOR avx+unroll version 
${fpath}/TEST './Table/dna3_40_4_16_new.table' './Input/dna/small.in' 40 4 8 8 16 2 8 >> RECORD

# FOR Unrolling and AVX speculation version
${fpath}/TEST './Table/dna3_40_4_16_new.table' './Input/dna/small.in' 40 4 8 8 16 2 10 >> RECORD

done

