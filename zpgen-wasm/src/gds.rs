// GDS file format export with streaming support
// Port from patternFileUtils.cpp

use std::io::{self, Write};

#[derive(Debug)]
pub enum GdsError {
    IoError(io::Error),
    InvalidData(String),
}

impl From<io::Error> for GdsError {
    fn from(err: io::Error) -> Self {
        GdsError::IoError(err)
    }
}

/// Trait for streaming GDS file output
pub trait GdsWriter {
    /// Write raw bytes to output
    fn write_chunk(&mut self, data: &[u8]) -> Result<(), GdsError>;

    /// Flush any buffered data
    fn flush(&mut self) -> Result<(), GdsError>;
}

/// GDS file writer that accumulates to a Vec<u8>
/// This is used for non-WASM environments and testing
pub struct VecGdsWriter {
    buffer: Vec<u8>,
}

impl VecGdsWriter {
    pub fn new() -> Self {
        Self {
            buffer: Vec::new(),
        }
    }

    pub fn into_inner(self) -> Vec<u8> {
        self.buffer
    }
}

impl Default for VecGdsWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl GdsWriter for VecGdsWriter {
    fn write_chunk(&mut self, data: &[u8]) -> Result<(), GdsError> {
        self.buffer.extend_from_slice(data);
        Ok(())
    }

    fn flush(&mut self) -> Result<(), GdsError> {
        Ok(())
    }
}

/// GDS header information
pub struct GdsHeader {
    pub layer_number: i32,
}

/// Initialize GDS file header
///
/// Port from patternFileUtils.cpp lines 123-152
pub fn init_gds<W: GdsWriter>(writer: &mut W, layer_number: i32) -> Result<(), GdsError> {
    #[rustfmt::skip]
    let gds_preamble: [u8; 102] = [
        0, 6, 0, 2, 0, 7, 0, 28, 1, 2, 230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
        230, 43, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 10, 2, 6, 110, 111, 110, 97,
        109, 101, 0, 20, 3, 5, 61, 104, 219, 139, 172, 113, 12, 180, 56, 109,
        243, 127, 103, 94, 246, 236, 0, 28, 5, 2, 0, 114, 0, 4, 0, 17, 0, 13,
        0, 22, 0, 56, 0, 114, 0, 4, 0, 17, 0, 13, 0, 22, 0, 56,
        0, 10, 6, 6, 110, 111, 110, 97, 109, 101
    ];

    writer.write_chunk(&gds_preamble)?;

    // Store layer number for polygon exports (would be in context in real implementation)
    Ok(())
}

/// Encode a 32-bit coordinate value into 4 bytes (big-endian)
///
/// Port from patternFileUtils.cpp lines 342-348
#[inline]
fn encode32(coord: i32) -> [u8; 4] {
    [
        ((coord >> 24) & 0xFF) as u8,
        ((coord >> 16) & 0xFF) as u8,
        ((coord >> 8) & 0xFF) as u8,
        (coord & 0xFF) as u8,
    ]
}

/// Export a polygon to GDS format
///
/// Port from patternFileUtils.cpp lines 362-395
pub fn export_polygon<W: GdsWriter>(
    writer: &mut W,
    coords: &[f64],
    layer_number: i32,
    num_vertices: usize,
) -> Result<(), GdsError> {
    // Each vertex has x,y coords = 2 values, 4 bytes each
    let num_coords = num_vertices * 2;
    let coord_buffer_size = num_coords * 4;
    let record_length = 4 + coord_buffer_size;

    // Polygon preamble: {0, 4, 8, 0, 0, 6, 13, 2, 0, [layer], 0, 6, 14, 2, 0, 0}
    let poly_pre: [u8; 16] = [
        0,
        4,
        8,
        0,
        0,
        6,
        13,
        2,
        0,
        layer_number as u8,
        0,
        6,
        14,
        2,
        0,
        0,
    ];

    // Polygon postamble: {0, 4, 17, 0}
    let poly_post: [u8; 4] = [0, 4, 17, 0];

    // Dynamic polygon format header with correct record length
    let poly_form: [u8; 4] = [
        ((record_length >> 8) & 0xFF) as u8,
        (record_length & 0xFF) as u8,
        16, // XY record type
        3,  // Four-byte signed integer
    ];

    // Write polygon structure
    writer.write_chunk(&poly_pre)?;
    writer.write_chunk(&poly_form)?;

    // Write coordinates
    for i in 0..num_coords {
        let coord_value = coords[i] as i32;
        let encoded = encode32(coord_value);
        writer.write_chunk(&encoded)?;
    }

    writer.write_chunk(&poly_post)?;

    Ok(())
}

/// Finalize GDS file
///
/// Port from patternFileUtils.cpp lines 154-158
pub fn render_gds<W: GdsWriter>(writer: &mut W) -> Result<(), GdsError> {
    let gds_postamble: [u8; 8] = [0, 4, 7, 0, 0, 4, 4, 0];

    writer.write_chunk(&gds_postamble)?;
    writer.flush()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode32() {
        let encoded = encode32(0x12345678);
        assert_eq!(encoded, [0x12, 0x34, 0x56, 0x78]);
    }

    #[test]
    fn test_gds_header() {
        let mut writer = VecGdsWriter::new();
        init_gds(&mut writer, 0).unwrap();

        let data = writer.into_inner();
        assert_eq!(data.len(), 102);
        assert_eq!(data[0], 0);
        assert_eq!(data[1], 6);
    }

    #[test]
    fn test_export_polygon() {
        let mut writer = VecGdsWriter::new();
        init_gds(&mut writer, 0).unwrap();

        // Simple square: 4 vertices + closing vertex = 5 vertices
        let coords = vec![
            0.0, 0.0, // v0
            100.0, 0.0, // v1
            100.0, 100.0, // v2
            0.0, 100.0, // v3
            0.0, 0.0, // v0 (close)
        ];

        export_polygon(&mut writer, &coords, 0, 5).unwrap();
        render_gds(&mut writer).unwrap();

        let data = writer.into_inner();
        assert!(data.len() > 102); // Should have header + polygon + footer
    }
}
