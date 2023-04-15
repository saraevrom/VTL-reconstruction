import numba as nb

# @nb.njit()
# def binsearch(array):
#     n = array.shape[0]
#     start = array[0]
#     end = array[-1]
#     if start>=0 and end>=0:
#         if start<=end:
#             return 0
#         else:
#             return n-1
#     elif start<=0 and end<=0:
#         if start>=end:
#             return 0
#         else:
#             return n-1
#     else:
#         # Binary search itself
#         start_i = 0
#         end_i = n-1
#         middle_i = (start_i+end_i)//2
#         while start_i!=middle_i:
#             if array[middle_i] == 0:
#                 return middle_i
#             if array[start_i]*array[middle_i]<0:
#                 end_i = middle_i
#             elif array[end_i]*array[middle_i]<0:
#                 start_i = middle_i
#             else:
#                 raise RuntimeError("How did we get there?")
#             middle_i = (start_i + end_i) // 2
#         return middle_i


@nb.njit(cache=True)
def binsearch_tgt(array, target):
    n = array.shape[0]
    start = target - array[0]
    end = target - array[-1]
    if start>=0 and end>=0:
        if start<=end:
            return 0
        else:
            return n-1
    elif start<=0 and end<=0:
        if start>=end:
            return 0
        else:
            return n-1
    else:
        # Binary search itself
        start_i = 0
        end_i = n-1
        middle_i = (start_i+end_i)//2
        while start_i != middle_i:
            mid_item = target - array[middle_i]
            start_item = target - array[start_i]
            end_item = target - array[end_i]
            if mid_item == 0:
                return middle_i
            if start_item*mid_item<0:
                end_i = middle_i
            elif end_item*mid_item<0:
                start_i = middle_i
            else:
                raise RuntimeError("How did we get there?")
            middle_i = (start_i + end_i) // 2
        return middle_i


@nb.njit(cache=True)
def binsearch(array):
    return binsearch_tgt(array, 0.0)

