from functools import cmp_to_key
import time


# Implementing the ternary split quicksort algorithm
# @param suffixes: A list of suffixes to be sorted
# @param depth: The current depth of the recursion
# @return: A sorted list of suffixes
def ternary_split_quicksort(suffixes, depth=0):
    if len(suffixes) <= 1:
        return suffixes

    # Determine the pivot character at the current depth
    pivot = suffixes[len(suffixes) // 2][depth] if depth < len(suffixes[len(suffixes) // 2]) else ''

    less = []
    equal = []
    greater = []

    for suffix in suffixes:
        if depth < len(suffix):
            char = suffix[depth]
            if char < pivot:
                less.append(suffix)
            elif char > pivot:
                greater.append(suffix)
            else:
                equal.append(suffix)
        else:
            # If depth exceeds the length of the suffix, consider it as 'less'
            less.append(suffix)

    # Recursively sort the partitions
    less_sorted = ternary_split_quicksort(less, depth)
    equal_sorted = ternary_split_quicksort(equal, depth + 1)
    greater_sorted = ternary_split_quicksort(greater, depth)

    return less_sorted + equal_sorted + greater_sorted

def construct_suffix_array(text):
    # Generate all suffixes
    suffixes = [text[i:] for i in range(len(text))]
    
    # Sort the suffixes using ternary split quicksort
    sorted_suffixes = ternary_split_quicksort(suffixes)

    # Construct the suffix array
    suffix_array = [len(text) - len(suffix) for suffix in sorted_suffixes]

    return suffix_array


def construct_suffix_array_custom(text, sorting_algorithm="insertion"):
    # Generate all suffixes
    suffixes = [text[i:] for i in range(len(text))]

    
    # Sort suffixes using bucket sort
    sorted_suffixes = bucket_sort_suffixes(suffixes, sorting_algorithm)

    # Construct the suffix array
    suffix_array = [len(text) - len(suffix) for suffix in sorted_suffixes]

    return suffix_array



def insertion_sort_suffixes(arr):
    for i in range(1, len(arr)):
        key = arr[i]
        j = i - 1
        while j >= 0 and arr[j] > key:
            arr[j + 1] = arr[j]
            j -= 1
        arr[j + 1] = key
    return arr



def selection_sort_suffixes(arr):
    for i in range(len(arr)):
        min_idx = i
        for j in range(i + 1, len(arr)):
            if arr[j] < arr[min_idx]:
                min_idx = j
        arr[i], arr[min_idx] = arr[min_idx], arr[i]
    return arr

def bucket_sort_suffixes(suffixes, sorting_algorithm="insertion"):
    if len(suffixes) == 0:
        return suffixes
    
    # Determine the number of buckets
    num_buckets = len(suffixes)
    max_suffix = max(suffixes)
    min_suffix = min(suffixes)

    # Create buckets and distribute the elements
    buckets = [[] for _ in range(num_buckets)]
    for suffix in suffixes:
        index = int((ord(suffix[0]) - ord(min_suffix[0])) / (ord(max_suffix[0]) - ord(min_suffix[0]) + 1) * num_buckets)
        buckets[index].append(suffix)

    # Sort each bucket using the specified sorting algorithm
    if sorting_algorithm == "insertion":
        for i in range(num_buckets):
            buckets[i] = insertion_sort_suffixes(buckets[i])
    elif sorting_algorithm == "selection":
        for i in range(num_buckets):
            buckets[i] = selection_sort_suffixes(buckets[i])
    else:
        raise ValueError("Invalid sorting algorithm or not implemented yet")

    sorted_suffixes = []
    for bucket in buckets:
        sorted_suffixes.extend(bucket)

    return sorted_suffixes


## The rest of the code is from the second mandatory assignment
def cmp_suffix(i,j):
    if (i < j):
        # Suffix i is lexicographically less than suffix j
        return -1
    if (i == j):
        # Suffix i and j are equal
        return 0
    else:
        # Suffix i is lexicographically greater than suffix j
        return 1
    
def suffixes(s):
    arr = []
    for i in range(len(s)):
        suffix = s[i:]
        arr.append(suffix)
    return arr



# Testing 
text = "banana"
suffix_array_builtin = sorted(suffixes(text), key = cmp_to_key(cmp_suffix))
print("Suffix array using built in sorting function:", suffix_array_builtin)

suffix_array_quick_sort = construct_suffix_array(text)
print("Suffix array using ternary split quicksort:", suffix_array_quick_sort)

suffix_array_insertion = construct_suffix_array_custom(text, sorting_algorithm="insertion")
print("Suffix array using insertion sort:", suffix_array_insertion)

suffix_array_selection = construct_suffix_array_custom(text, sorting_algorithm="selection")
print("Suffix array using selection sort:", suffix_array_selection)


######################################
#### REAL TESTING FOR REAL MEN #######
#### With random genetics string #####
######################################
x = 'gactaagttaacacacgcatcggccgttccgacatccaactggttcctccatcggggccaatttccattctactgagggattccttgaggatattgcgagcgcgcgtggacattagcggtgtgtctgaaccacgtcttagacgtcaatatctgcgcaatcaagaaatgatcgtccaactgaggatgcaccatcttgtataacgcaatagacggactcatggtaagtgtcagcgtagtatcacagaactatcagagatctaaatcctcagttccttgtgagttggctaagaccagcgccggtactccggggggccgtatcttaggtgtaaccccgaagttcgccccagcttaagcagccacttagatgagggcctcgctgccagacgtcctcgctgggggccaatttataaccgactcagttacatccgcggggagttcgttgaaacaccggaggccgctgggtagtctttgtcgtacgtactttcggagcttccattcgggtgctcgacacagcgcgcttaaactcgctgttcataacacctggatatacgtgcgtcgagctgctattcttccatatacgatagggctagacacatcgcttaagtggccatcccttagacttctttacctgtgcgttagcttcattcggtttgcaaacgcaggcagtttgagccgtgaactggacgagaacgcgtggtgtctggcatattctttccggttcggaatagagggtcttttgatcgattcctgggccgtgtggctcgaaactttctatggagccgccagcagtaaagatgcatttacttgtgatatgggcaatccgcttctcgttggatccaggcaataaaaaagttccctcccgattcagtattcggactgagataggccatgaggtgatgttcagacttctaggatcgcggacgctgacccgcatattcgaccacggagacggcgtaatccgtgaatgtccttatttcgcactgttacacgtgctgagaggtgttgccgccaaaaaccgctgggttcaacagccggcaccggtttcgacatcgtttttcagcggcctcactaattattgtacccgagcaggaacagattgcagcgggctgccacgcgaaccaacctggagggtggtgtgggcggttcaatgacgtagcctctatcaagactgaaaggaaagtatgtacgatctctagccagacacgccgttatacttgagttaaggcttcaccctctcgagaacccgcgatgattggtatcgcgccgttcaattctctataaagattcttattcacagccccaggccgtggtacgtcaagcggatgcggaataaaaccacacaaaggattttccggggtccgagttcgtatcacgtatggtagaggttagaaatattttgtgaattgcaattctggcgataatcgttgcagtctgatggcggcgtgattagggtcggcatccctcgcagaaatgggagggcgcctcctgaccgtctaacggttatctgaaacggatgttagcgcgaagcatctaatatatcttcgggtacgactggaggtcaagtgggccgactggcggctttaactgcttaggcatgtatagccattgccaagaaatgccaggctgccaggggaggacctcagagagcggaccgtttgttatcccagagaggaggggacagggacacttctccaaattcgtccggttcataaagaccacttttatacatatggctctatcggttaggtgagggaggctgatactgtttgcgatcgtacatctgacctgtgagttcccgttaggaattatccgcttgtagattttccggacagacacaaaatgtataatgtccgtactccatagtaaaaccctatcattacttcatggagccgggcccgagtttccctgcgagaagccttgacctcctgagattagcgatccatcctaatatgagatcccaagcctgacatatggaccagattcggtcactgactctgatgttcacatctaatgagtaggttacgccccgagacacaactctaacaggacccacctaaacgtcgtacggttggacctggttcatacagctagtgtctctcaagggctcatggacaaccccgcaggataaagtatgtgtccaatttatcgcacggtcgttcttttcaataaagcccattagatctgaattgtttctgactttatcatagtgcaagtggatgatgacgagtggggccctcgacgataatagatgcctgctactcatacgccggtggaggcgagtgaactgacccggtcgcgtataggctgcggtagaggtttcttgaatgtagccacaaatgcaatgacatcatctctctagatcacgtcagctgaatacggaagtcgatgaaccgccaaagattcctacttcaaatcgtagcctatttttagttccttgggtaagacagtcgtgaccgaaagcaggtatatgcatctggtatccatttgccacgcacctagcaccccgctatagtttaatcgctcaattaccctcctagaatcgaagtctgagaaggctaggttagtaaattggccgctgtaggcggtgcatgcgcacaggcatttgcctatcgctaactgaggtccagacagttatacagcagactcataataaccgcacccaccgcagtagccctatctatctaccccattcttttcgttgaaggttccgattctgccttgtagccggagacgacatggcttcctgacttgctataataacgtcgccgagggacatgacccgggactctgagggctcaatacttgtgagtcgtccgccagtacggacatgaggattctacaaatcctgataaaagatgtacgcgtcgatgccccccgctaagcgcatagtgatctgtactctaaccaagtaaactgtgcagcatcgtgatcaccccggaaactgcgttccgcgtatgtgcctaaaggcttcgcagagatatcgccatcgctttccggcgggggtctccggccttaacggggcttcagaaaaatgcaatttctctatcgttgaacgccgggtgttccaaagaattcggaatcaaagcctgcgaggacaagggtttacctataaagccaccaggactaatcgattggcgtgacgctgaaaagttcgaaggtgcaggcttcgtctcgcattaggtgtcatcaagtggaattaaaaggcggggcgggcgtcaaagtcgttgggctcctacatactgcaagtattacggtggttggggtaagcctcactatacggcctgtcggccgcctcttactaggcattccttccttaaagtacaccttattcgagcacacacacccgtactacacccgaggtctgtctgtacagacatggtgccatttgctggtgtaatacggcccactacatagggacaaggcatcctccgcgagtctaccaaatactgcgatttctatcctatggaatttcggacggtcgatgggcggaacaagaaggaagatgcctacacaacatcaataggtccaattcaattgctctctcatcgctatgagtgtggttaaaggtctcctatttagaattagaaggcacgggcagtatttccctctgcgatttcgtcacgcttaaagtcatcccgatatgttgtacgagtctagagtaatccgcgcattgtaatgcttacacttgagtcagaagggaggtgtactaggtcatatttctccacgattgttccccgtgaaggcgtagttaccaaatatgtaatagataatgtcatagggggactgtcctctaagtaccctcttctaacgcataaccatgagtgattctggccacaagatatctctgcgaggtcagtcagtaattccatcaattacggtgtgcaacaagagacagacactaatagtttaggtctactatgtcgggagggatttctttggagcctgtgggcaattcgtgaaatggacctaactaagcttacttgcagaggcattgttcctgagctccccttaattaattcaaaccagagatgacagttgtacttaacgtagaatcagcagttgaaggatacaatcttttggaagtagttactagttgccacgatgcagtctgaagtttattgagcgaacggacactagcggatatgtataagacgtccaccgtcgcagcttcgcacttattaattaaggaggcgcttgacattttaacgagtattctctaccctaaaactctgttagaatcgcgagcctaagcaaaaagcggcgggatctatcgaaggtaatacactcataatgaagtagtccgggtgcgggggacacaggccacatcctaccgagccatacctgtgcatcatagtagacagagtgactgataagtgatcagtcattagtttgccactacttgtcttactcgacccaggagcgaatactcggagctcggtagcgctctttcatgattttggacggcacatggcgagttaagctaggatcatgcttcgtaccacccccgctatgagatctcgtcaattgggactactgtctttagtcagctcggataactaaaagtggggacaatgccggatgtggcttaaggctacgatcgtatccagggatgtcaggtctcctgtcataatgcgatcagagtacagtctcatgctaacgagggacgggaaccaagttcaatgctggcgtattggcctcacccccagaggtgtcccaggatgtgtatctataaatttcagactcactcgaaagaacttgatactcttgccgtggcggtcaagctttgcgatcctactcccctgtgataatgcctacgtccgcctaagctgctggtgcgatggactcccgacggctcggtactactgatagcgctatatcgacagtgttaataagtctgagccccttcatagataattccaggtcaccagcgagaacgatagcggccgagtccccatgacatcaacgtgggctcctctgtgcctgcgaagtaagcctagtatgctgtcggtgtattacgaaccgatttggagctcagttaagttccctagtgtgatttcacttagcaccttcgaacgctctccacaccaattatggtagttgtttcagtagtttcactagattccatgagggatctgcggtcagtagcaggatcggaaatcacccgtcatgtttaagcctataaggacgttcagcggtttgggatgcttaatgaggttatacgacacatgattcctggtctgccgtagtgtacaacataagtaacttgttacgctctccaccttttagaatattgtcaggccccaacgggataacggaacgtggaacagatccgtagcaaagaggcaatcttccgcgaatgtgaagcgaacccactcgcatgatacaattctatgggcataccattcccaaacattaagtataacgctggatttaggggaacaaggtctgacactgagtgaacagcgccagacaattaaaacggatgaaacgacaattaactagtgctaacagcagatatagattagggatactctttttgtgccacccatcgtaatatccaaaactctccggttcaagccctttcaactactttccgctcccgccgtgccaaaacaacagggcggtgattacgaacatgattatgttccgacgcttcgaagaactcaaagttagctgtagctagttcgcccagtcccgaaacagtcggcccttgagctcgagcgaccagtatgtggcaaaattccgaactatagtactactctgtcgttaacattttcctataaatcgtgaagcttagctccccttttagcggccgtaacgagaataatagttcaagaaataccgcattcgcgtgtcatctggttcgcgccctcgcgtgttcgacagccagacacaaatactccagtagggagcgacgtcagagtcggatcccagagtcatctcgtgaagctgtcagtgctggacaggtacagtcaagacaattaaatcgtctagagcctacccagattggagtatactggtgtgccaataagggcttgcctacaccctcgcgagaacagagtgtgtaacttcgaaatgccgtcgagcccttataaaaggattgtctaggattttccgagtctcgatctggccgattgcaaagtgacgcctgtatgaggcgcgttatacactaacatgcctgagataaagaatacttcgggatgccgcgcgtctgatggtaaggtagcgatcgggaggcctaccggaccgaatgaaccgtatgaataactaaatccccatctgattcaacggtaaactcattatgcagccgagggggcgtcacattttttcgataaactcgtacgtggctccgtggaagctgtctcaggatcctggtagaatctcccccctcctcttaccctaaaccgggagatggcccctataataagggactatgtaactatcaagagtttgaagtggtacccaaactcatgcccatcccttgtgtacaactaagggtattacacatgcatatgccactagaaatacgtcgatgattcgtcgagtagagtttcttgccttgggaacttgcgcttggattgtctgtcacgcgggacgtgtctttgcccatgacaggtctgcgtgacttagatacgatcttctagcggttgcacgtgtaccgtaatttgatatacgtggccacatacgttcgtaactcatgtctggacgggattataaaagggtaggacactatctcgatccacgggtagtgcagcaacacgatggtcattaaggatggcccggacgacgttgcaaaggcgccacaccacgctcgctcttaactagagatcggcagctgagaccgcgtacaggatattattaagcctgactgagacaaagcagcgatggacctagcggtgatttgacccaagcatatagcgcatgtttacagacgccaccacccatactattacaaggtctcctcactttgaaacccaacacgccctttctagacatccagagtcgattccaaacggtaaaaagggaggggaaccttaacagctacgaccaactgcgatgagcgactccagaaagcatcagttccccacgatattgaatttgagtatatgttcccatggacaggttacaatagaccttactaacaagcgctttatatcacctcagcccatttcttctctccctaccgataggctcccgtctatggtccgccagtgccttgttcgaagtaggtaccgctatttccgcccggcgagtattactcaggaacctaccttacagatgtcacggtgccatcaagccacagcatggcccaagagcttttgtcgtcaggccgtgactttacgctactccctttatgagcgcatatactcgagaatgtcgcatccgctcgtctggataccgaggttccgtgattctctttattgcaagggcagagaatgtggcatgatcgcccctgatcgacatagcctccccattcggtacttcaccgctctacgtaatctcgtgcatgactcaggcggaaattgacgagcacaaccagagtgtatcctggagattcacggactatgcgcggctgttataggatgagtacgcaccttcaaccagtcgaatgttaaaccggttgcgacagccctgtctagcgggttataaaagcttgttggatacggaacccccctcaggtatattaccccccaaagatacgcagtttcgaattacggtcggtgatctaattccgcaccccggtgcgtagggcgcagcggcagatcccaaatctcttcgtgaggccacagattataggctaacttgctctataatcctgcccgctgtctctagaccaaacgatgacttatttagcgagtccgcactagcagggccctcgctttccgcaaaaccacggcgcacgtaaagggtagtcgggagccaggcgagtctctaccacatctgaatggagctaggggtagccgacgttgtcgactttcacaaaatcgcacgtttctcttgacagtatttggtctcgcttgcggcttggatattaccctgcatctatgagacctgcgaataacgcgggaaatgtgcccacaacgcccaacgttctctacaccctctcggcaaccgcccgactgggagtcatcagcttgtagagtcgatgagttatcggactacataggccttggcacattttaaaattgagggttagttcaagtcccgcaaccatcttagaccacgtcactattatgcgtagcgagcatacgcggcaatgatcccgtcagtaacagtcgtgaataacgcagcggaatcacgagcggcatcgccatctatgtactatgagtggaaaacggtgcggcacgtatttccctaggtatatacttagtcaaccgctatggtttacgtaaatggggttgcgcaggacatcgtaagatgactaaccatcccaccggcatatggagcctttgtttaaagcatgttggttgtcgtatggcagaactctccccgccgtagcccctgatatgcaacgcatgattgactacaacataagtcagaccacctgacatcaaacagcgcccaacgattcttgtgaacaatcccctgttcccacaaagggaagcaccgcgcagggatgcataaccaatacgctgtccggatagacttctacgtttgctacggggccttttaggttatagccattattactgtgatgtgcggaaggtgcagtagaagttggtcatagacttagattggagagtacatcgaccgaaagacttcaaccagtgtgcaggtagttttagtcgacgtaagttagttcagagacttcaatcaccttaacatagaccatataagagattaggttcttcaaattgccttcggtgattgaggcgtcacggcctgatcgatgtgacccgttgcttaaccgccgacactcaattcaccactcgtctcggactcttctgaaacgcgtacctgcgacggacccatgctcatgaaacacacaaatccggttatgtttttctccttcctaactagtatcacctgctaaaactaaacgtatcattccgccgtgcgattaagaagctgaaatgaagtggcgtaacgagtgaataacttctttaaagtaaggagcctgggtcgagatcaggaccacgatgcgggtccagacgcagccaggatcttgcgacgtataacgcctttaagtctaggctgttggactccacttccatcctgtagacgcgcattcttggattgaagctggaggtcataactcgtattgaaggctcttgccccctcttgccgcacagaagaaagggccctagtatttggattagaattgttgcgttaatacactcgtgaaaactgcgatagctcctctattagtcggcggcgcgcgaaaaaaactcgttacggcttcctcagagacttcggacacacctcagccaccggtttgttcgtgcacgtgacacagggtcccgctgggttcgagatcctattacatacataacttattatgccctaccttaggccgctcagtaggctgggcttaccgcaataacccttgttagtgctctaggaagggtatatagcgagaactgtagactcctgctcggtattgtaagtcgggaaacttgtcgaaggcagtatccggattttagtgggcctcggcccgtgtttcccatttcaactaaagacggagccggagctaacgtagccgcgtgcatcctttattacgtgggagcatcccgatgcgatctatagtagtcggatagaatcgtaagattctagtcaccatgcgtggagattaggatagcccacaaccccggcactgattcaaatgtcaaaagatttcgtttccattgtaggcttggattcacgaccgccctaaacatttagcattcagggccggggcgccgacaacaatcgttagtcagaagggaacctgggcctacgactccaggttcgagaccgggtccttgggaagagcgttttcctgtggaaaatagaaatgatc'
start_time_builtin = time.time()
suffix_array_builtin = sorted(suffixes(x), key = cmp_to_key(cmp_suffix))
end_time_builtin = time.time()
running_time_builtin = end_time_builtin - start_time_builtin
print(f'Running time with string is {round(running_time_builtin, 4)} seconds using built in sorting function')

start_time_quick_sort = time.time()
suffix_array_quick_sort = construct_suffix_array(x)
end_time_quick_sort = time.time()
running_time_quick_sort = end_time_quick_sort - start_time_quick_sort
print(f'Running time with string is {round(running_time_quick_sort, 4)} seconds using ternary split quicksort')

start_time_insertion = time.time()
suffix_array_insertion = construct_suffix_array_custom(x, sorting_algorithm="insertion")
end_time_insertion = time.time()
running_time_insertion = end_time_insertion - start_time_insertion
print(f'Running time with string is {round(running_time_insertion, 4)} seconds using insertion sort')

start_time_selection = time.time()
suffix_array_selection = construct_suffix_array_custom(x, sorting_algorithm="selection")
end_time_selection = time.time()
running_time_selection = end_time_selection - start_time_selection
print(f'Running time with string is {round(running_time_selection, 4)} seconds using selection sort')



######################################
#### REAL TESTING FOR REAL MEN #######
#### With a string of a's      #######
######################################
x = 'a'*10000
start_time_builtin = time.time()
suffix_array_builtin = sorted(suffixes(x), key = cmp_to_key(cmp_suffix))
end_time_builtin = time.time()
running_time_builtin = end_time_builtin - start_time_builtin
print(f'Running time with string is {round(running_time_builtin, 4)} seconds using built in sorting function')

start_time_quick_sort = time.time()
suffix_array_quick_sort = construct_suffix_array(x)
end_time_quick_sort = time.time()
running_time_quick_sort = end_time_quick_sort - start_time_quick_sort
print(f'Running time with string is {round(running_time_quick_sort, 4)} seconds using ternary split quicksort')

start_time_insertion = time.time()
suffix_array_insertion = construct_suffix_array_custom(x, sorting_algorithm="insertion")
end_time_insertion = time.time()
running_time_insertion = end_time_insertion - start_time_insertion
print(f'Running time with string is {round(running_time_insertion, 4)} seconds using insertion sort')

start_time_selection = time.time()
suffix_array_selection = construct_suffix_array_custom(x, sorting_algorithm="selection")
end_time_selection = time.time()
running_time_selection = end_time_selection - start_time_selection
print(f'Running time with string is {round(running_time_selection, 4)} seconds using selection sort')
