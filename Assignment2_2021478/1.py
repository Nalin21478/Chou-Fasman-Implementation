seq="SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"
alpha_dictionary = {'A': 1.45, 'R': 0.79, 'N': 0.73, 'D': 0.98, 'C': 0.77, 'E': 1.53, 'Q': 1.17, 'G': 0.53, 'H': 1.24, 'I': 1.00, 'L': 1.34, 'K': 1.07, 'M': 1.20, 'F': 1.12, 'P': 0.59, 'S': 0.79, 'T': 0.82, 'W': 1.14, 'Y': 0.61, 'V': 1.14}
beta_dictionary = {'A': 0.97, 'R': 0.90, 'N': 0.65, 'D': 0.80, 'C': 1.30, 'E': 0.26, 'Q': 1.23, 'G': 0.81, 'H': 0.71, 'I': 1.60, 'L': 1.22, 'K': 0.74, 'M': 1.67, 'F': 1.28, 'P': 0.62, 'S': 0.72, 'T': 1.20, 'W': 1.19, 'Y': 1.29, 'V': 1.65}
# Initialize lists for H Positions, S positions, residue and final sequence
residue = []
final_seq = []
H_arr = []
S_arr = []
print("For the given sequence:\n"+ seq+"\nThe output of CHou-Fasman is as follows:\n")

'''The following function identifies, wehre exactly the H/S should come and then it returns an array that stores the updated positions'''

seq_len=len(seq)
def find_position(seq, compare_dict):
    arr = []
    for i in range(seq_len-5):
        seq_change = seq[i:i+6]
        count = 0
        for x in seq_change:
            if compare_dict.get(x) >= 1:
                count+=1
        if count >= 4:
            for j in range(6):
                if not(i+j in arr):
                    arr.append(i+j)
            exit_cond = 1000000
            exitss = i+6
            while exit_cond >= 4:
                if exitss < seq_len:
                    extended_seq=seq[exitss-3:exitss+1]
                    exit_cond = 0
                    for y in extended_seq:
                        exit_cond +=compare_dict.get(y)
                    if exit_cond >= 4:
                        if not(exitss in arr):
                            arr.append(exitss)
                else:
                    break
                exitss +=1
            exit_cond =(pow(10,6))
            exitss = i-1
            while exit_cond >= 4:
                if exitss >= 0:
                    extended_seq=seq[exitss:exitss+4]
                    exit_cond = 0
                    for y in extended_seq:
                        exit_cond = exit_cond +compare_dict.get(y)
                    if exit_cond >= 4:
                        if not(exitss in arr):
                            arr.append(exitss)
                else:
                    break
                exitss -=1
    return arr

H_arr=find_position(seq,alpha_dictionary)#For the H positions
S_arr=find_position(seq,beta_dictionary)#for the S positions


# Set all elements of final sequence to 0
for l in range(seq_len):
    final_seq.append(0)

# Update final sequence based on H_arr and S_arr
for i in range(seq_len):
    if (i in H_arr) and (i in S_arr):
        residue.append(i)
    elif (i in H_arr) and not (i in S_arr):
        final_seq[i] = 'H'
    elif not(i in H_arr) and (i in S_arr):
        final_seq[i] = 'S'
    else:
        final_seq[i] = '*'

# Initialize variables
i = 0
counterr = 0
a1 = 0
b1 = 0

# Update final sequence based on residue
while i < len(residue):
    vkk = i
    residualll = residue[i]
    i += 1
    if i >= len(residue):
        break
    while residue[i] == residualll + 1:
        residualll += 1
        i += 1
        if i >= len(residue):
            break
    counterr += 1
    extended_seq = seq[residue[vkk]:residue[i-1]+1]
    
    for x in extended_seq:
        a1 += alpha_dictionary.get(x)
        b1 += beta_dictionary.get(x)
    if a1 < b1:
        for j in range(i-vkk):
            final_seq[residue[vkk]+j] = 'S'
       
    else:
        for j in range(i-vkk):
            final_seq[residue[vkk]+j] = 'H'


# Print helix positions
H_pos_string=''
for i in range(seq_len):
    if(i not in H_arr):
        H_pos_string+='*'
    else:
        H_pos_string+='H'
print("The helix positions are as follows:\n"+H_pos_string+"\n")

# Print beta positions
S_pos_string=''
for j in range(seq_len):
    if(j not in S_arr):
        S_pos_string+='*'
    else:
        S_pos_string+='H'
print("The beta positions are as follows:\n"+S_pos_string+"\n")           

# Print final H and S positions
HNS_pos_string=''
for i in final_seq:
    HNS_pos_string+=i

print("The final H and S positions are as follows:\n"+HNS_pos_string)


            
        
    