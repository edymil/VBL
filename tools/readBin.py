#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 13:47:43 2017

@author: thierry
"""

import binascii
import struct

def text_to_bits(text, encoding='utf-8', errors='surrogatepass'):
    bits = bin(int(binascii.hexlify(text.encode(encoding, errors)), 16))[2:]
    return bits.zfill(8 * ((len(bits) + 7) // 8))

def text_from_bits(bits, encoding='utf-8', errors='surrogatepass'):
    n = int(bits, 2)
    return int2bytes(n).decode(encoding, errors)

def int2bytes(i):
    hex_string = '%x' % i
    n = len(hex_string)
    return binascii.unhexlify(hex_string.zfill(n + (n & 1)))
  
  
from functools import partial

with open("0.bin", mode='rb') as file: # b is important -> binary
    fileContent = file.read()
    print(struct.unpack("iiiii", fileContent[:20]))
    print(struct.unpack("i" * ((len(fileContent) -24) // 4), fileContent[20:-4]))
    print(struct.unpack("i", fileContent[-4:]))
#with open("0.bin", "rb") as f:
#    byte = f.read(1)
#    while byte != "":
#        # Do stuff with byte.
#        byte = f.read(1)
#        a=binascii.b2a_hqx(byte)
#with open("0.bin", "rb") as f:
#    while True:
#        byte = f.read(1)
#        if not byte:
#            break
#        text_from_bits(ord(byte))
#with open('0.bin', 'rb') as file:
#    for byte in iter(partial(file.read, 1), b''):
#        text_from_bits(byte)