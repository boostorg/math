from sympy import symbols, simplify, ccode
from sympy.tensor.array import derive_by_array
import sys

def generate_cpp_tensor(expr, vars, order):
    derivatives = expr
    for _ in range(order):
        derivatives = derive_by_array(derivatives, vars)

    # Flatten tensor for variable assignment
    def flatten(tensor):
        if not hasattr(tensor, '__iter__') or isinstance(tensor, (str, bytes)):
            return [tensor]
        if hasattr(tensor, 'tolist'):
            tensor = tensor.tolist()
        flat_list = []
        for e in tensor:
            flat_list.extend(flatten(e))
        return flat_list

    flat_derivs = flatten(derivatives)

    # Generate multi-indices for all tensor elements
    def generate_indices(d, order):
        if order == 0:
            return [[]]
        else:
            smaller = generate_indices(d, order - 1)
            return [s + [i] for s in smaller for i in range(d)]

    all_indices = generate_indices(len(vars), order)

    # Generate variable names like f_x, f_xy, f_xyz, etc.
    var_names = []
    for idx in all_indices:
        suffix = ''.join(str(vars[i]) for i in idx)
        var_names.append(f"f_{suffix}")

    # Assign each derivative to a separate variable
    code_lines = []
    for var_name, expr in zip(var_names, flat_derivs):
        simplified = simplify(expr)
        c_expr = ccode(simplified)
        code_lines.append(f"    T {var_name} = static_cast<T>({c_expr});")

    # Now build nested vector initialization recursively matching the shape of derivatives
    def build_nested_vector(tensor, level=0, index_start=0):
        if level == order:
            # At the scalar level, return the variable name at flat index
            return var_names[index_start], 1
        else:
            # tensor is iterable, build vector of sub-elements
            if hasattr(tensor, 'tolist'):
                tensor = tensor.tolist()
            elems = []
            count = 0
            for t in tensor:
                s, c = build_nested_vector(t, level+1, index_start + count)
                elems.append(s)
                count += c
            return '{ ' + ', '.join(elems) + ' }', count

    nested_vector_str, _ = build_nested_vector(derivatives)

    # Compose return type string depending on order
    def vector_type(level):
        if level == 0:
            return "T"
        else:
            return f"std::vector<{vector_type(level-1)}>"

    return_type = vector_type(order)

    # Compose final C++ function body with separate variable assignments and nested return vector
    body = '\n'.join(code_lines) + f'\n    return {nested_vector_str};'

    return return_type, body

if __name__ == "__main__":
    x, y, z = symbols('x y z')
    vars = [x, y]
    f = x*x*x*x*y*y*y*y 
    order = int(sys.argv[1]) if len(sys.argv) > 1 else 2

    print(f"// Order-{order} derivative of f(x, y) = {f}")

    ret_type, func_body = generate_cpp_tensor(f, vars, order)

    print(f"template<typename T>")
    print(f"{ret_type} gf_a(T x, T y) {{")
    print(func_body)
    print("}")
