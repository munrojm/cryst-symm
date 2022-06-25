use crate::structure::Structure;

#[derive(Debug)]
pub struct Reducer {
    pub dtol: f32,
}

impl Reducer {
    pub fn find_primitive_cell(structure: Structure) -> Structure {
        return structure;
    }
}
