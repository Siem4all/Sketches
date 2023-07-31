import math
import random


class CntrMaster:
    def __init__(self, cntrSize, hyperSize):
        self.cntrSize = cntrSize
        self.hyperSize = hyperSize
        self.exponentMax = 2 ** self.hyperSize - 1
        self.exponent_range = [[sum([2 ** i for i in range(self.exponentMax - j + 1, self.exponentMax)]),
                                sum([2 ** i for i in range(self.exponentMax - j, self.exponentMax)])] for j in
                               range(self.exponentMax)]

    def calculateExponentSizeUsingExponentValue(self, exponent_value):
        """
        Returns the size of the exponent field given the value of the hyper-exponent field.
        """
        for j in self.exponent_range:
            if j[0] <= exponent_value < j[1]:
                return self.exponentMax - self.exponent_range.index(j)

    def exponentValue(self, exponent_bits, exponent_size):
        """
        Returns the size of the exponent field given the value of the hyper-exponent field.
        """
        return self.maxBiasValue() - sum([(1 + exponent_bits[i]) for i in range(exponent_size)])

    def CalculateMantissaSizeUsingExponentValue(self, exponent_value):
        """
        Returns the size of the mantissa field given the size the counter,
        the size of the hyper-exponent field, and the size of the exponent field.
        """
        return self.cntrSize - self.hyperSize - self.calculateExponentSizeUsingExponentValue(exponent_value)

    def mantissaValue(self, mantissa_bits, mantissa_size):
        """
        Returns the size of the mantissa value given the mantissa bits and mantissa size.
        """
        return sum([mantissa_bits[i] * (2 ** i) for i in range(mantissa_size)])

    def offsetValue(self, exponentValue):
        """
        Returns the offset value given the exponent value and mantissa size.
        """
        return sum([2 ** (i + self.CalculateMantissaSizeUsingExponentValue(i)) for i in range(exponentValue)])

    def maxBiasValue(self):
        """
        Returns the maximal bias value using exponent maximum value.
        """
        return sum([2 ** i for i in range(1, self.exponentMax)])

    def counterMaxValue(self):
        """
        Returns the max value represented by the counter.
        """
        return sum(
            [2 ** (i + int(self.CalculateMantissaSizeUsingExponentValue(i))) for i in range(self.maxBiasValue())]) + (
                       2 ** (self.cntrSize - self.hyperSize) - 1) * (2 ** int(self.maxBiasValue()))

    def counterValue(self, exponentValue, mantissaValue):
        """
        Returns the value represented by the counter given the mantissa and exponent fields.
        """
        return self.offsetValue(exponentValue) + mantissaValue * 2 * exponentValue

    # Define the IncCounter algorithm
    def IncCounter(self, T, Ec):
        En = min([E for E in range(Ec, self.maxBiasValue() + 1) if self.offsetValue(E) <= T])
        Mn = (T - self.offsetValue(En)) // (2 ** En)
        if self.counterValue(En, Mn) < self.counterMaxValue():
            if Mn < (2 ** self.CalculateMantissaSizeUsingExponentValue(En)) - 1:
                p = (T - self.counterValue(En, Mn)) / (self.counterValue(En, Mn + 1) - self.counterValue(En, Mn))
                if random.random() < p:
                    Mn += 1
            else:
                p = (T - self.counterValue(En, Mn)) / (self.counterValue(En + 1, Mn) - self.counterValue(En, Mn))
                if random.random() < p:
                    En += 1
                    Mn = 0
        return En, Mn


# Define the increment_counter function
def main():
    f2p = CntrMaster(4, 2)
    counter = 0
    Ec = 0
    increment = 100
    target = counter + increment
    En, Mn = f2p.IncCounter(target, Ec)

    Tl = f2p.counterValue(En, Mn)

    En, Mn = f2p.IncCounter(target, En)

    Th = f2p.counterValue(En, Mn)

    print(Tl, Th)


# Test the increment_counter function
if __name__ == '__main__':
    main()
