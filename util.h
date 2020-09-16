#define RED           "\033[31;1m"
#define GREEN         "\033[32;1m"
#define YELLOW        "\033[33;1m"
#define BLUE          "\033[34;1m"
#define PURPLE        "\033[35;1m"
#define CYAN          "\033[36;1m"
#define WHITE         "\033[37;1m"
#define URED          "\033[31;1;4m"
#define UGREEN        "\033[32;1;4m"
#define UYELLOW       "\033[33;1;4m"
#define UBLUE         "\033[34;1;4m"
#define UPURPLE       "\033[35;1;4m"
#define UCYAN         "\033[36;1;4m"
#define UWHITE        "\033[37;1;4m"
#define DEFAULTCOLOR "\033[0m"

#ifndef _PARALEL_
#define COLOR(color, arg) color << arg << DEFAULTCOLOR
#else
#define COLOR(color, arg) arg
#endif