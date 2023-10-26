// #ifndef DEBUGOUT_H
// #define DEBUGOUT_H

// #ifdef LOGGING
// uint32_t logging = 0;

// void debug(clause *clause, uint32_t loglevel) {
//   if (loglevel > logging)
//     return;
//   std::cout << "c DEBUG CLAUSE:  ";
//   for (auto lit : *clause) {
//     std::cout << lit << ", ";
//   }
//   std::cout << std::endl;
// }

// #define dout0 std::cout << "c DEBUG (" << __LINE__ << "): "
// #define dout1                                                                  \
//   if (logging > 0)                                                             \
//   std::cout << "c DEBUG (" << __LINE__ << "): "
// #define dout2                                                                  \
//   if (logging > 1)                                                             \
//   std::cout << "c DEBUG (" << __LINE__ << "): "
// #define dout3                                                                  \
//   if (logging > 2)                                                             \
//   std::cout << "c DEBUG (" << __LINE__ << "): "

// #else
// #define dout0 0 && std::cout
// #define dout1 0 && std::cout
// #define dout2 0 && std::cout
// #define dout3 0 && std::cout
// #define debug(...)                                                             \
//   do {                                                                         \
//   } while (0)
// #endif

// #endif // DEBUGOUT_H