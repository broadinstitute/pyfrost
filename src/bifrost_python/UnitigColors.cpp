#include "UnitigColors.h"

namespace pyfrost {

UnitigColorIterator::UnitigColorIterator()
{
}

UnitigColorIterator::UnitigColorIterator(UnitigColors::const_iterator const& to_wrap) :
    iter(to_wrap)
{
    if(!iter.isInvalid()) {
        current = iter.getColorID();
    }
}

UnitigColorIterator::value_type UnitigColorIterator::operator*() const {
    return current;
}

UnitigColorIterator::pointer UnitigColorIterator::operator->() const {
    return &current;
}

void UnitigColorIterator::operator++() {
    iter.nextColor();
    if(!iter.isInvalid()) {
        current = iter.getColorID();
    }
}

UnitigColorIterator UnitigColorIterator::operator++(int) {
    UnitigColorIterator tmp(*this);
    operator++();

    return tmp;
}

bool UnitigColorIterator::operator==(UnitigColorIterator const& o) {
    return iter == o.iter;
}

bool UnitigColorIterator::operator!=(const UnitigColorIterator &o) {
    return !operator==(o);
}


void define_UnitigColors(py::module &m) {
    auto py_UnitigColors = py::class_<UnitigColorsProxy>(m, "UnitigColors")
        .def("__getitem__", &UnitigColorsProxy::getColorsAtPos, py::is_operator())
        .def("__contains__", &UnitigColorsProxy::contains)
        .def("__len__", &UnitigColorsProxy::size)
        .def("__iter__", [](UnitigColorsProxy const &self) {
            return py::make_iterator(self.begin(), self.end());
        })

        .def_static("_from_iterable", [] (py::iterable const& iterable) {
            return py::set(iterable);
        })

        .def("num_kmers_with_color", &UnitigColorsProxy::numKmersWithColor);

    auto Set = py::module::import("collections.abc").attr("Set");
    py_UnitigColors.attr("__bases__") = py::make_tuple(Set).attr("__add__")(py_UnitigColors.attr("__bases__"));
}

}
