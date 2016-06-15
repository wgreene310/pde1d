#include <stdexcept>
#include <string>

class PDE1dException : public std::exception {
public:
  PDE1dException(const char *id, const char *msg) :
    id(id), msg(msg) {}
  const char *getId() const { return id.data();  }
#if __cplusplus <= 199711L
  virtual const char *what() const { return msg.data(); }
#else
  virtual const char *what() const noexcept{ return msg.data(); }
#endif
private:
  std::string id, msg;
};