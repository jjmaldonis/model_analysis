from collections import Counter
import pyvoro

def compute_indexes(coords, xsize=1000, ysize=1000, zsize=1000, index=None):
    result = pyvoro.compute_voronoi(coords,
        [[-xsize/2., xsize/2], [-ysize/2., ysize/2], [-zsize/2., zsize/2]], # limits
        len(coords), # block size
    )
    if index:
        return convert_result_to_index(result[index])
    else:
        for i, r in enumerate(result):
            result[i] = convert_result_to_index(r)
        return result

def convert_result_to_index(result):
    index = [len(f['vertices']) for f in result['faces']]
    index = Counter(index)
    index = tuple(index[i] for i in range(3,13))
    return index




def main():
    coords = [[-0.499999823344, 0,  0.809016708539],
         [0.499999823344,  0,  -0.809016708539],
         [0.0, -0.809016708539, 0.499999823344],
         [0.0, 0.809016708539,  -0.499999823344],
         [0.0, 0.809016708539,  0.499999823344],
         [0.0, -0.809016708539, -0.499999823344],
         [0.809016708539,  0.499999823344,  -0],
         [0.809016708539,  -0.499999823344, 0],
         [0.499999823344,  -0,  0.809016708539],
         [-0.809016708539, -0.499999823344, 0.0],
         [-0.809016708539, 0.499999823344,  -0],
         [-0.499999823344, 0,  -0.809016708539],
         [0,  0,  0]
    ]
    index = compute_indexes(coords, index=-1)
    print(index)



if __name__ == '__main__':
    main()
