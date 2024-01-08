#ifndef ELEMENT_TYPE_H
#define ELEMENT_TYPE_H

#include <vector>
#include <flowFormat.hpp>
#include <string>
#include <map>


class elementTypeFormat
{
public:
    geom_int id;
    std::string name;

    geom_int nNodes;
    std::vector<geom_int> nodesOrder;
    std::vector<std::vector<geom_int>>  nodesOrderPlanes;

    std::map<int, std::string> elementNameFromGmshID = {
        {2, "triangle"}, {3, "quad"}, {4, "tetra"}, {5, "hex"}, {6, "prism"}, {7, "pyramid"},
    };

    elementTypeFormat() {};

    elementTypeFormat( std::string &_name, geom_int _nNodes,// std::vector<geom_int> &_nodesOrder,
                 std::vector<std::vector<geom_int>> &_nodesOrderPlanes ) 
    {
        //this->id   = _id;
        this->name = _name;
        this->nNodes = _nNodes;
//        this->nodesOrder = _nodesOrder;
        this->nodesOrderPlanes = _nodesOrderPlanes;
    };
};

class elementTypeMap
{
public:
    std::map<std::string, elementTypeFormat> mapElementFromName;
    std::map<geom_int, elementTypeFormat>    mapElementFromGmshID;

    elementTypeMap()
    {
        // Quad
        std::string _name3 = "triangle";
        geom_int _nNodes3  = 3;
        std::vector<std::vector<geom_int>> _nodesOrderFaces3{ {0, 1},
                                                              {1, 2}, 
                                                              {2, 0}, }; 

        elementTypeFormat triangle = elementTypeFormat( _name3, _nNodes3, _nodesOrderFaces3);


        // Quad
        std::string _name4 = "quad";
        geom_int _nNodes4  = 4;
        std::vector<std::vector<geom_int>> _nodesOrderFaces4{ {0, 1},
                                                              {1, 2}, 
                                                              {2, 3}, 
                                                              {3, 0}, }; 

        elementTypeFormat quad = elementTypeFormat( _name4, _nNodes4, _nodesOrderFaces4);

        // Hexahedral
        std::string _name5 = "hex";
        geom_int _nNodes5  = 8;
        std::vector<std::vector<geom_int>> _nodesOrderFaces5{ {1, 2, 6, 5},   // u+
                                                              {0, 4, 7, 3},   // u-
                                                              {2, 3, 7, 6},   // v+
                                                              {0, 1, 5, 4},   // v-
                                                              {4, 5, 6, 7},   // w+
                                                              {0, 3, 2, 1},}; // w-

        elementTypeFormat hex = elementTypeFormat(_name5, _nNodes5, _nodesOrderFaces5);

        // tetra
        std::string _name6 = "tetra";
        geom_int _nNodes6  = 4;
        std::vector<std::vector<geom_int>> _nodesOrderFaces6{ {0, 3, 2},   
                                                              {1, 2, 3},   
                                                              {0, 2, 1},   
                                                              {0, 1, 3},};   

        elementTypeFormat tetra = elementTypeFormat(_name6, _nNodes6, _nodesOrderFaces6);

        // prism
        std::string _name7 = "prism";
        geom_int _nNodes7  = 6;
        std::vector<std::vector<geom_int>> _nodesOrderFaces7{ {3, 4, 5},   
                                                              {0, 2, 1},   
                                                              {0, 1, 4, 3},   
                                                              {1, 2, 5, 4},
                                                              {0, 3, 5, 2},
                                                              };   

        elementTypeFormat prism = elementTypeFormat(_name7, _nNodes7, _nodesOrderFaces7);

        // pyramid
        std::string _name8 = "pyramid";
        geom_int _nNodes8  = 5;
        std::vector<std::vector<geom_int>> _nodesOrderFaces8{ {0, 1, 4},   
                                                              {0, 4, 3},   
                                                              {3, 4, 2},   
                                                              {2, 4, 1},   
                                                              {0, 3, 2, 1},   
                                                              };   

        elementTypeFormat pyramid = elementTypeFormat(_name8, _nNodes8, _nodesOrderFaces8);



        mapElementFromName.insert(std::make_pair("triangle", triangle));
        mapElementFromName.insert(std::make_pair("quad", quad));
        mapElementFromName.insert(std::make_pair("hex" , hex));
        mapElementFromName.insert(std::make_pair("tetra" , tetra));
        mapElementFromName.insert(std::make_pair("prism" , prism));
        mapElementFromName.insert(std::make_pair("pyramid" , pyramid));

        mapElementFromGmshID.insert(std::make_pair(2, triangle));
        mapElementFromGmshID.insert(std::make_pair(3, quad));
        mapElementFromGmshID.insert(std::make_pair(4, tetra));
        mapElementFromGmshID.insert(std::make_pair(5, hex));
        mapElementFromGmshID.insert(std::make_pair(6, prism));
        mapElementFromGmshID.insert(std::make_pair(7, pyramid));
    }

};

/* Copied from https://gmsh.info/doc/texinfo/gmsh.html

Line:                  Line3:           Line4:

      v
      ^
      |
      |
0-----+-----1 --> u    0----2----1      0---2---3---1

Triangle:               Triangle6:          Triangle9/10:

v
^
|
2                       2                    2
|`\                     |`\                  | \
|  `\                   |  `\                7   6
|    `\                 5    `4              |     \
|      `\               |      `\            8  (9)  5
|        `\             |        `\          |         \
0----------1--> u       0-----3----1         0---3---4---1


Triangle12/15:

 v
 ^
 |
 2
 | \
 9   8
 |     \
10 (14)  7
 |         \
11 (12) (13) 6
 |             \
 0---3---4---5---1--> u


Quadrangle:            Quadrangle8:            Quadrangle9:

      v
      ^
      |
3-----------2          3-----6-----2           3-----6-----2
|     |     |          |           |           |           |
|     |     |          |           |           |           |
|     +---- | --> u    7           5           7     8     5
|           |          |           |           |           |
|           |          |           |           |           |
0-----------1          0-----4-----1           0-----4-----1

Tetrahedron:                          Tetrahedron10:

                   v
                 .
               ,/
              /
           2                                     2
         ,/|`\                                 ,/|`\
       ,/  |  `\                             ,/  |  `\
     ,/    '.   `\                         ,6    '.   `5
   ,/       |     `\                     ,/       8     `\
 ,/         |       `\                 ,/         |       `\
0-----------'.--------1 --> u         0--------4--'.--------1
 `\.         |      ,/                 `\.         |      ,/
    `\.      |    ,/                      `\.      |    ,9
       `\.   '. ,/                           `7.   '. ,/
          `\. |/                                `\. |/
             `3                                    `3
                `\.
                   ` w


Hexahedron:             Hexahedron20:          Hexahedron27:

       v
3----------2            3----13----2           3----13----2
|\     ^   |\           |\         |\          |\         |\
| \    |   | \          | 15       | 14        |15    24  | 14
|  \   |   |  \         9  \       11 \        9  \ 20    11 \
|   7------+---6        |   7----19+---6       |   7----19+---6
|   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
 \  |    \  \  |         \  17      \  18       \ 17    25 \  18
  \ |     \  \ |         10 |        12|        10 |  21    12|
   \|      w  \|           \|         \|          \|         \|
    4----------5            4----16----5           4----16----5


Prism:                      Prism15:               Prism18:

           w
           ^
           |
           3                       3                      3
         ,/|`\                   ,/|`\                  ,/|`\
       ,/  |  `\               12  |  13              12  |  13
     ,/    |    `\           ,/    |    `\          ,/    |    `\
    4------+------5         4------14-----5        4------14-----5
    |      |      |         |      8      |        |      8      |
    |    ,/|`\    |         |      |      |        |    ,/|`\    |
    |  ,/  |  `\  |         |      |      |        |  15  |  16  |
    |,/    |    `\|         |      |      |        |,/    |    `\|
   ,|      |      |\        10     |      11       10-----17-----11
 ,/ |      0      | `\      |      0      |        |      0      |
u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    |
    |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  |
    |,/         `\|         |,/         `\|        |,/         `\|
    1-------------2         1------9------2        1------9------2


Pyramid:                     Pyramid13:

               4                            4
             ,/|\                         ,/|\
           ,/ .'|\                      ,/ .'|\
         ,/   | | \                   ,/   | | \
       ,/    .' | `.                ,/    .' | `.
     ,/      |  '.  \             ,7      |  12  \
   ,/       .' w |   \          ,/       .'   |   \
 ,/         |  ^ |    \       ,/         9    |    11
0----------.'--|-3    `.     0--------6-.'----3    `.
 `\        |   |  `\    \      `\        |      `\    \
   `\     .'   +----`\ - \ -> v  `5     .'        10   \
     `\   |    `\     `\  \        `\   |           `\  \
       `\.'      `\     `\`          `\.'             `\`
          1----------------2            1--------8-------2
                    `\
                       u

Pyramid14:

               4
             ,/|\
           ,/ .'|\
         ,/   | | \
       ,/    .' | `.
     ,7      |  12  \
   ,/       .'   |   \
 ,/         9    |    11
0--------6-.'----3    `.
  `\        |      `\    \
    `5     .' 13     10   \
      `\   |           `\  \
        `\.'             `\`
           1--------8-------2
                    `\
                       u

*/

#endif