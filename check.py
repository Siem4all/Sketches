num = 255
n_bits = 8
cntrMaxVal  = (1 << n_bits) - 1  # the max counter value which is 2**cntrSize
bits = format(num, f'0{n_bits}b')
print(bits, cntrMaxVal, 2**n_bits)
