#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/**
 * Library-level functions.
 * You should use them in the main sections.
 */

uint64_t convertToUint64 (double number) {
    return *((uint64_t *)(&number));
}

bool getBit (const uint64_t number, const uint8_t index) {
    uint64_t tmp = number;
    tmp >>= index;
    return tmp;
}

uint64_t invertBit (const uint64_t number, const uint8_t index) {
    return getBit (number, index)
        ? number ^ (1 << index)
            : number | (1 << index);
}


/**
 * Checkers here:
 */

bool checkForPlusZero (uint64_t number) {
    return number == 0x0000000000000000;
}

bool checkForMinusZero (uint64_t number) {
    return number == 0x8000000000000000;
}

bool checkForPlusInf (uint64_t number) {
    return number == 0x7ff0000000000000;
}

bool checkForMinusInf (uint64_t number) {
    return number == 0xfff0000000000000;
}

bool not_null_or_inf (uint64_t number) {
    return !checkForPlusZero(number) && !checkForMinusZero(number)
            && !checkForPlusInf(number) && !checkForMinusInf(number);
}

bool checkForSignalingNan (uint64_t number) {
    return not_null_or_inf(number)
           && ((number & 0x7ff8000000000000) == 0x7ff0000000000000);
}

bool checkForQuietNan (uint64_t number) {
    return not_null_or_inf(number)
            && ((number & 0x7ff8000000000000) == 0x7ff8000000000000);
}

bool checkForPlusNormal (uint64_t number) {
    return not_null_or_inf(number)
            && !checkForQuietNan(number) && !checkForSignalingNan(number)
            && ((number & 0xfff8000000000000) < 0x7fff000000000000)
            && ((number & 0xfff8000000000000) > 0x0008000000000000);
}

bool checkForMinusNormal (uint64_t number) {
    return not_null_or_inf(number)
            && !checkForQuietNan(number) && !checkForSignalingNan(number)
            && ((number & 0xfff8000000000000) < 0xfff0000000000000)
            && ((number & 0xfff8000000000000) > 0x0008000000000000);
}

bool checkForPlusDenormal (uint64_t number) {
    return not_null_or_inf(number)
           && !checkForQuietNan(number) && !checkForSignalingNan(number)
           && (!(number & 0xfff0000000000000) < 0x0010000000000000);
}

bool checkForMinusDenormal (uint64_t number) {
    return not_null_or_inf(number)
           && !checkForQuietNan(number) && !checkForSignalingNan(number)
           && ((number & 0xfff0000000000000) == 0x800000000000000);
}

void classify (double number) {
    if (checkForPlusZero(convertToUint64(number))) {
        printf("Plus zero\n");
    }

    else if (checkForMinusZero(convertToUint64(number))) {
        printf("Minus zero\n");
    }

    else if (checkForPlusInf(convertToUint64(number))) {
        printf("Plus inf\n");
    }

    else if (checkForMinusInf(convertToUint64(number))) {
        printf("Minus inf\n");
    }

    else if (checkForPlusNormal(convertToUint64(number))) {
        printf("Plus normal\n");
    }

    else if (checkForMinusNormal(convertToUint64(number))) {
        printf("Minus normal\n");
    }

    else if (checkForPlusDenormal(convertToUint64(number))) {
        printf("Plus Denormal\n");
    }

    else if (checkForMinusDenormal(convertToUint64(number))) {
        printf("Minus Denormal\n");
    }

    else if (checkForSignalingNan(convertToUint64(number))) {
        printf("Signailing NaN\n");
    }

    else if (checkForQuietNan(convertToUint64(number))) {
        printf("Quiet NaN\n");
    }

    else {
        printf("Error.\n");
    }
}

int main() {

    double input;
    scanf("%lf", &input);
    classify(input);

    return 0;
}
