#include "NodeDataDict.h"
#include "UnitigColors.h"

namespace pyfrost {

void define_NodeDataDict(py::module& m) {

    auto handle = py::class_<NodeDataDict>(m, "NodeDataDict")
        // Dict like interface
        .def("__getitem__", &NodeDataDict::getData)
        .def("__setitem__", &NodeDataDict::setData)
        .def("__delitem__", &NodeDataDict::delData)
        .def("__contains__", &NodeDataDict::contains, py::is_operator())
        .def("__len__", &NodeDataDict::size)
        .def("__iter__", [](NodeDataDict const& self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        .def("__str__", &NodeDataDict::mappedSequence)

        .def("twin", [](NodeDataDict const& self) {
            return self.unitigTail().twin();
        })
        .def("rep", py::overload_cast<>(&NodeDataDict::unitigRepresentative, py::const_),
            "Return the head k-mer of this unitig in forward strand")
        .def("full_node", [](NodeDataDict const& self) {
            return NodeDataDict(self.mappingToFullUnitig());
        });

    // Hack to let our NodeDataDict inherit from the collections.abc.MutableMapping Python mixin
    auto MutableMapping = py::module::import("collections.abc").attr("MutableMapping");
    handle.attr("__bases__") = py::make_tuple(MutableMapping).attr("__add__")(handle.attr("__bases__"));
}

const std::unordered_map<std::string, UnitigMetaKeys> NodeDataDict::hardcoded_keys = {
    {"length",          UnitigMetaKeys::LENGTH},
    {"pos",             UnitigMetaKeys::POS},
    {"strand",          UnitigMetaKeys::STRAND},
    {"head",            UnitigMetaKeys::HEAD},
    {"tail",            UnitigMetaKeys::TAIL},
    {"mapped_sequence", UnitigMetaKeys::MAPPED_SEQUENCE},
    {"unitig_sequence", UnitigMetaKeys::UNITIG_SEQUENCE},
    {"unitig_length",   UnitigMetaKeys::UNITIG_LENGTH},
    {"is_full_mapping", UnitigMetaKeys::IS_FULL_MAPPING},
    {"colors",          UnitigMetaKeys::COLORS},
};

bool NodeDataDict::contains(const std::string &key) const {
    if(getMetaKey(key) == UnitigMetaKeys::NONE) {
        return getDataDict().contains(key.c_str());
    }

    return true;
}

py::object NodeDataDict::getData(std::string const& key) const {
    switch(getMetaKey(key)) {
        case UnitigMetaKeys::LENGTH:
            return py::cast(unitig.len);

        case UnitigMetaKeys::POS:
            return py::cast(unitig.strand
                ? unitig.dist
                : unitig.dist > 0 ? (unitig.size - unitig.getGraph()->getK() - unitig.dist) : 0);

        case UnitigMetaKeys::STRAND:
            return py::cast(unitig.strand ? Strand::FORWARD : Strand::REVERSE);

        case UnitigMetaKeys::HEAD:
            return py::cast(unitigHead());

        case UnitigMetaKeys::TAIL:
            return py::cast(unitigTail());

        case UnitigMetaKeys::MAPPED_SEQUENCE:
            return py::cast(mappedSequence());

        case UnitigMetaKeys::UNITIG_SEQUENCE:
            return py::cast(unitigSequence());

        case UnitigMetaKeys::UNITIG_LENGTH:
            return py::cast(unitigLength());

        case UnitigMetaKeys::IS_FULL_MAPPING:
            return py::cast(unitig.isFullMapping());

        case UnitigMetaKeys::COLORS:
        {
            if(!unitig.strand) {
                auto um = unitig;

                // Transform to forward orientation if necessary
                um.strand = true;
                um.dist = 0;

                return py::cast(UnitigColorsProxy(um, true));
            } else {
                return py::cast(UnitigColorsProxy(unitig));
            }
        }

        default:
        {
            auto& data_dict = getDataDict();

            if(data_dict.contains(key)) {
                return data_dict[key.c_str()];
            } else {
                throw py::key_error("Key '" + key + "' not found in unitig metadata.");
            }
        }
    }
}

void NodeDataDict::setData(std::string const& key, py::object const& value) {
    if(getMetaKey(key) != UnitigMetaKeys::NONE) {
        throw py::key_error("Key '" + key + "' is read only unitig metadata.");
    }

    auto& data_dict = getDataDict();
    data_dict[key.c_str()] = value;
}

void NodeDataDict::delData(const py::str &key) {
    if(getMetaKey(key) != UnitigMetaKeys::NONE) {
        throw py::key_error("Key '" + key.cast<std::string>() + "' is read only unitig metadata.");
    }

    auto& data_dict = getDataDict();
    PyDict_DelItem(data_dict.ptr(), key.ptr());
}

size_t NodeDataDict::size() const {
    return hardcoded_keys.size() + getDataDict().size();
}

std::string NodeDataDict::mappedSequence() const {
    return unitig.mappedSequenceToString();
}

size_t NodeDataDict::mappedSequenceLength() const {
    return unitig.getGraph()->getK() + unitig.len - 1;
}

std::string NodeDataDict::unitigSequence() const {
    return unitig.strand ? unitig.referenceUnitigToString() : reverse_complement(unitig.referenceUnitigToString());
}

size_t NodeDataDict::unitigLength() const {
    return unitig.size - unitig.getGraph()->getK() + 1;
}

Kmer NodeDataDict::unitigHead() const {
    return unitig.strand ? unitig.getUnitigHead() : unitig.getUnitigTail().twin();
}

Kmer NodeDataDict::unitigTail() const {
    return unitig.strand ? unitig.getUnitigTail() : unitig.getUnitigHead().twin();
}

Kmer NodeDataDict::unitigRepresentative() const {
    return unitig.getUnitigHead();
}

UnitigDataKeyIterator NodeDataDict::begin() const {
    auto& data_dict = getDataDict();

    return {hardcoded_keys.cbegin(), hardcoded_keys.cend(),
        data_dict.begin()};
}

UnitigDataKeyIterator NodeDataDict::end() const {
    auto& data_dict = getDataDict();
    return {hardcoded_keys.cend(), hardcoded_keys.cend(),
            data_dict.end()};
}

UnitigMetaKeys NodeDataDict::getMetaKey(std::string const &key) const {
    auto key_enum = hardcoded_keys.find(key);
    if(key_enum != hardcoded_keys.end()) {
        return key_enum->second;
    }

    return UnitigMetaKeys::NONE;
}

}
