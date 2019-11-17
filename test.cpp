#include <memory>
#include <iostream>

using namespace std;

int main() {
    int *a = new int;
    *a = 0;
    unique_ptr<int> u(a);
    *u = 1;

    cout << *a << " " << *u << endl;
    cout << a << " " << &(*u) << endl;
}