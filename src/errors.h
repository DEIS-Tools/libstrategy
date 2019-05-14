/*
 *  Copyright Peter G. Jensen, all rights reserved.
 */

/* 
 * File:   errors.h
 * Author: Peter G. Jensen <root@petergjoel.dk>
 *
 * Created on May 10, 2019, 9:46 AM
 */

#ifndef ERRORS_H
#define ERRORS_H

struct base_error : public std::exception {
    std::string _message;

    explicit base_error(std::string m)
    : _message(std::move(m)) {
    }

    const char *what() const noexcept override {
        return _message.c_str();
    }

    virtual void print(std::ostream &os) const {
        os << what() << std::endl;
    }

    friend std::ostream &operator<<(std::ostream &os, const base_error &el) {
        el.print(os);
        return os;
    }
};

#endif /* ERRORS_H */

