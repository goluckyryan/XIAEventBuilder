/************************************
Clover      :   0 -  99
BGO         : 100 - 199 
GAGG        : 200 - 299
ZERO DEGREE : 300 - 399

 * *********************************/
#ifndef MAPPING
#define MAPPING

//==================== mapping

#define NCLOVER        10
#define NCRYSTAL       NCLOVER*4
#define NBGO           NCLOVER
#define NGAGG          26
#define NZERO          2

//                    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
int mapping[128] ={   0,   1,   2,   3, 100,   4,   5,   6,   7, 101,   8,   9,  10,  11, 102,  -1,  //mod-0
                     12,  13,  14,  15, 103,  16,  17,  18,  19, 104,  20,  21,  22,  23, 105,  -1,  //mod-1
                     24,  25,  26,  27, 106,  28,  29,  30,  31, 107,  32,  33,  34,  35, 108,  -1,  //mod-2
                     36,  37,  38,  39, 109,  40,  41,  42,  43, 110, 300, 301, 200, 201, 202, 203,  
                    204, 205, 206, 207, 208, 209, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259,   
                    210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225,  // Ring 4A
                    260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275,  // Ring 4B
                     -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1};  

//200- 209 GAGG 2A
//210- 225 GAGG 4A

//250- 259 GAGG 2B
//260- 275 GAGG 4B

                 
#endif
