#ifndef INLINES_H
#define INLINES_H

inline int calculate_box_number(const int i, const int j, const int k, const ivec3 N_boxes)
{
//    if (i > N_boxes(0) || j > N_boxes(1) || k > N_boxes(2))
//    {
//        cout << "!!! calculate_box_number: Box index larger than number of boxes. Aborting!" << endl;
//        abort();
//    }
    return i + j*N_boxes(0) + k*N_boxes(0)*N_boxes(1);
}

inline int calculate_box_number(const ivec box_index, const ivec3 N_boxes)
{
//    if ((box_index(0) > N_boxes(0)) != 0 || (box_index(1) > N_boxes(1)) != 0 || (box_index(2) > N_boxes(2)) != 0)
//    {
//        cout << "!!! calculate_box_number: Box index larger than number of boxes. Aborting!" << endl;
//        abort();
//    }
    return box_index(0) + box_index(1)*N_boxes(0) + box_index(2)*N_boxes(0)*N_boxes(1);
}

inline ivec3 calculate_box_index(const int box_number, const ivec3 N_boxes)
{
    ivec3 box_index;
    box_index << box_number%N_boxes(0)
              << (box_number/N_boxes(0))%N_boxes(1)
              << ((box_number/N_boxes(0))/N_boxes(1))%N_boxes(2);
             // << (box_number/N_boxes(0))/N_boxes(1); // same as the above line
    return box_index;
}

#endif // INLINES_H
