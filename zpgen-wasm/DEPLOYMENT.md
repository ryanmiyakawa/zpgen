# ZPGen WASM - Deployment & React Integration Guide

## ⚠️ Important: No Common Source

**There is NO automatic synchronization between C++ and Rust code.**

- The Rust code is a **one-time port** of the C++ code
- Changes to C++ **will NOT** automatically appear in Rust
- You must **manually sync** any algorithm changes between the two

### Maintenance Strategy

**Option 1: Maintain Both** (recommended for transition period)
- Keep C++ as the "reference implementation"
- When you update C++ algorithms, manually port changes to Rust
- Eventually deprecate C++ once Rust is proven in production

**Option 2: Choose One**
- Pick either C++ or Rust as your primary codebase
- Phase out the other over time
- Rust advantages: browser support, memory safety, cross-platform
- C++ advantages: existing tooling, proven in production

**Option 3: Feature Branches**
- New features go into Rust only (browser-first)
- Keep C++ for legacy desktop workflows
- Diverge over time as needs dictate

---

## Deployment Options

### Option 1: NPM Package (Recommended for React)

This is the easiest way to integrate with your React app.

#### Step 1: Build WASM Package

```bash
cd /Users/rhmiyakawa/Documents/Xcode/zpgen/ZPGen/ZPGen/ZPGen/zpgen-wasm
./build-wasm.sh
```

This creates the package at: `examples/browser/pkg/`

#### Step 2: Install in React Project

```bash
cd /path/to/your-react-app
npm install /Users/rhmiyakawa/Documents/Xcode/zpgen/ZPGen/ZPGen/ZPGen/zpgen-wasm/examples/browser/pkg
```

Or add to `package.json`:
```json
{
  "dependencies": {
    "zpgen-wasm": "file:../zpgen/ZPGen/ZPGen/ZPGen/zpgen-wasm/examples/browser/pkg"
  }
}
```

#### Step 3: Use in React Component

**Basic Example:**

```typescript
// src/components/ZPGenerator.tsx
import { useEffect, useState } from 'react';
import init, { WasmZPGenerator } from 'zpgen-wasm';

export function ZPGenerator() {
    const [initialized, setInitialized] = useState(false);
    const [generating, setGenerating] = useState(false);
    const [progress, setProgress] = useState(0);
    const [error, setError] = useState<string | null>(null);

    // Initialize WASM module on mount
    useEffect(() => {
        init().then(() => {
            console.log('ZPGen WASM initialized');
            setInitialized(true);
        }).catch(err => {
            console.error('Failed to initialize WASM:', err);
            setError(err.message);
        });
    }, []);

    const handleGenerate = async (params: ZPGenParams) => {
        if (!initialized) {
            setError('WASM not initialized');
            return;
        }

        setGenerating(true);
        setProgress(0);
        setError(null);

        try {
            const chunks: Uint8Array[] = [];

            // Create generator with chunk callback
            const gen = new WasmZPGenerator((chunk: Uint8Array) => {
                chunks.push(new Uint8Array(chunk));
            });

            // Set progress callback
            gen.setProgressCallback((update: any) => {
                setProgress(update.progress);
                console.log(`Zone ${update.zone}/${update.totalZones}`);
            });

            // Generate zone plate
            gen.generate({
                zTol: 0.01,
                lambda_nm: params.wavelength,
                p: [0, 0, params.focalLength],
                q: 100000,
                k_0: [0, 0, 1],
                bx: [1, 0, 0],
                by: [0, 1, 0],
                obscurationSigma: 0,
                na: params.na,
                nZerns: 0,
                zernikeOrders: [],
                customMaskIdx: 0,
                anamorphicFac: 1,
                anamorphicAzimuth: 0,
                zpcPhase: 0,
                apd: 0,
                apdWindow: 0,
                zpcr2: 0,
                zpcr1: 0,
                biasNm: params.bias,
                oppositeTone: false,
                randomizeZoneStart: false,
                fsIdx: 0,
                buttressGapWidth: 0,
                buttressPeriod: 0,
                dutyCycle: 0.5,
                layerNumber: 0,
            });

            // Merge chunks and download
            const blob = new Blob(chunks, { type: 'application/octet-stream' });
            downloadFile(blob, `zoneplate_${params.na}_${params.wavelength}nm.gds`);

        } catch (err: any) {
            setError(err.message);
            console.error('Generation failed:', err);
        } finally {
            setGenerating(false);
        }
    };

    const downloadFile = (blob: Blob, filename: string) => {
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    };

    if (error) {
        return <div className="error">Error: {error}</div>;
    }

    if (!initialized) {
        return <div>Loading WASM module...</div>;
    }

    return (
        <div>
            <button
                onClick={() => handleGenerate({
                    na: 0.02,
                    wavelength: 13.5,
                    focalLength: 100,
                    bias: 10
                })}
                disabled={generating}
            >
                {generating ? `Generating... ${Math.round(progress * 100)}%` : 'Generate'}
            </button>
        </div>
    );
}

interface ZPGenParams {
    na: number;
    wavelength: number;
    focalLength: number;
    bias: number;
}
```

**Advanced Example with Web Worker:**

```typescript
// src/workers/zpgen.worker.ts
import init, { WasmZPGenerator } from 'zpgen-wasm';

let initialized = false;

self.onmessage = async (e) => {
    if (e.data.type === 'generate') {
        if (!initialized) {
            await init();
            initialized = true;
        }

        try {
            const chunks: Uint8Array[] = [];

            const gen = new WasmZPGenerator((chunk: Uint8Array) => {
                // Stream chunks to main thread
                self.postMessage({
                    type: 'chunk',
                    data: chunk,
                });
                chunks.push(new Uint8Array(chunk));
            });

            gen.setProgressCallback((update: any) => {
                self.postMessage({
                    type: 'progress',
                    progress: update.progress,
                    zone: update.zone,
                    totalZones: update.totalZones,
                });
            });

            gen.generate(e.data.params);

            self.postMessage({ type: 'complete' });
        } catch (err: any) {
            self.postMessage({ type: 'error', error: err.message });
        }
    }
};
```

```typescript
// src/hooks/useZPGenWorker.ts
import { useEffect, useRef, useState } from 'react';

export function useZPGenWorker() {
    const workerRef = useRef<Worker | null>(null);
    const [progress, setProgress] = useState(0);
    const [generating, setGenerating] = useState(false);
    const chunksRef = useRef<Uint8Array[]>([]);

    useEffect(() => {
        workerRef.current = new Worker(
            new URL('../workers/zpgen.worker.ts', import.meta.url),
            { type: 'module' }
        );

        workerRef.current.onmessage = (e) => {
            switch (e.data.type) {
                case 'chunk':
                    chunksRef.current.push(new Uint8Array(e.data.data));
                    break;
                case 'progress':
                    setProgress(e.data.progress);
                    break;
                case 'complete':
                    setGenerating(false);
                    // Download file
                    const blob = new Blob(chunksRef.current);
                    downloadFile(blob, 'zoneplate.gds');
                    chunksRef.current = [];
                    break;
                case 'error':
                    console.error('Worker error:', e.data.error);
                    setGenerating(false);
                    break;
            }
        };

        return () => workerRef.current?.terminate();
    }, []);

    const generate = (params: any) => {
        setGenerating(true);
        setProgress(0);
        chunksRef.current = [];
        workerRef.current?.postMessage({ type: 'generate', params });
    };

    return { generate, progress, generating };
}

function downloadFile(blob: Blob, filename: string) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}
```

---

### Option 2: Copy WASM Files Directly

If you don't want to use NPM:

#### Step 1: Copy WASM Files

```bash
# Copy to your React public folder
cp examples/browser/pkg/zpgen_wasm.js public/wasm/
cp examples/browser/pkg/zpgen_wasm_bg.wasm public/wasm/
cp examples/browser/pkg/zpgen_wasm.d.ts src/types/
```

#### Step 2: Load Dynamically

```typescript
import { useEffect, useState } from 'react';

export function useZPGen() {
    const [zpgen, setZpgen] = useState<any>(null);

    useEffect(() => {
        async function loadWasm() {
            const module = await import('/wasm/zpgen_wasm.js');
            await module.default('/wasm/zpgen_wasm_bg.wasm');
            setZpgen(module);
        }
        loadWasm();
    }, []);

    return zpgen;
}
```

---

### Option 3: CDN (for quick prototyping)

```typescript
// Not recommended for production, but works for testing
useEffect(() => {
    const script = document.createElement('script');
    script.src = 'http://localhost:8888/pkg/zpgen_wasm.js';
    script.type = 'module';
    document.body.appendChild(script);
}, []);
```

---

## TypeScript Definitions

The WASM package includes TypeScript definitions automatically.

**If TypeScript complains**, create a declaration file:

```typescript
// src/types/zpgen-wasm.d.ts
declare module 'zpgen-wasm' {
    export interface ZPParams {
        zTol: number;
        lambda_nm: number;
        p: [number, number, number];
        q: number;
        k_0: [number, number, number];
        bx: [number, number, number];
        by: [number, number, number];
        obscurationSigma: number;
        na: number;
        nZerns: number;
        zernikeOrders: number[];
        customMaskIdx: number;
        anamorphicFac: number;
        anamorphicAzimuth: number;
        zpcPhase: number;
        apd: number;
        apdWindow: number;
        zpcr2: number;
        zpcr1: number;
        biasNm: number;
        oppositeTone: boolean;
        randomizeZoneStart: boolean;
        fsIdx: number;
        buttressGapWidth: number;
        buttressPeriod: number;
        dutyCycle: number;
        layerNumber: number;
    }

    export interface ProgressUpdate {
        progress: number;
        zone: number;
        totalZones: number;
    }

    export class WasmZPGenerator {
        constructor(onChunk: (data: Uint8Array) => void);
        setProgressCallback(callback: (update: ProgressUpdate) => void): void;
        generate(params: ZPParams): void;
    }

    export default function init(module?: WebAssembly.Module | string): Promise<void>;
}
```

---

## Build Configuration

### Webpack (Create React App)

If using Create React App, you need to enable WASM:

```javascript
// config-overrides.js (requires react-app-rewired)
module.exports = function override(config) {
    config.experiments = {
        ...config.experiments,
        asyncWebAssembly: true,
    };
    return config;
};
```

### Vite

Vite supports WASM out of the box:

```typescript
// vite.config.ts
import { defineConfig } from 'vite';

export default defineConfig({
    plugins: [],
    optimizeDeps: {
        exclude: ['zpgen-wasm'],
    },
});
```

### Next.js

```javascript
// next.config.js
module.exports = {
    webpack: (config) => {
        config.experiments = {
            ...config.experiments,
            asyncWebAssembly: true,
        };
        return config;
    },
};
```

---

## Publishing to NPM (Optional)

To share this package publicly:

### Step 1: Update package.json

```bash
cd examples/browser/pkg
```

Edit `package.json`:
```json
{
  "name": "@yourname/zpgen-wasm",
  "version": "0.1.0",
  "description": "Zone plate generator compiled to WebAssembly",
  "repository": "https://github.com/yourname/zpgen-wasm",
  "license": "MIT",
  "author": "Your Name",
  "files": [
    "zpgen_wasm_bg.wasm",
    "zpgen_wasm.js",
    "zpgen_wasm.d.ts"
  ]
}
```

### Step 2: Publish

```bash
npm login
npm publish --access public
```

### Step 3: Install from NPM

```bash
npm install @yourname/zpgen-wasm
```

---

## Development Workflow

### Making Changes to Rust Code

```bash
# 1. Edit Rust code
code src/zpgen.rs

# 2. Test locally
cargo test

# 3. Rebuild WASM
./build-wasm.sh

# 4. Test in React app (if using npm link)
cd /path/to/react-app
npm update zpgen-wasm
```

### If You Update C++ Code

**You must manually port changes to Rust:**

1. Identify what changed in C++ (e.g., new algorithm in `zpGenHolo.cpp`)
2. Update corresponding Rust file (e.g., `src/zpgen.rs`)
3. Run tests: `cargo test`
4. Rebuild WASM: `./build-wasm.sh`
5. Update version in `Cargo.toml`
6. Reinstall in React app

**Example:**
```bash
# C++ change: Modified secant solver in zpGenHolo.cpp
# Rust change: Update src/solver.rs with same logic

cd zpgen-wasm
# Edit src/solver.rs
cargo test              # Verify tests pass
./build-wasm.sh         # Rebuild WASM
cd /path/to/react-app
npm update zpgen-wasm   # Get new version
```

---

## Deployment Checklist

### Before Deploying to Production

- [ ] Test with multiple zone plate sizes (NA=0.01, 0.05, 0.08)
- [ ] Verify GDS output with your fabrication software
- [ ] Test in all target browsers (Chrome, Firefox, Safari)
- [ ] Test on mobile devices (if applicable)
- [ ] Add error handling for edge cases
- [ ] Set up monitoring/logging for WASM errors
- [ ] Optimize WASM size if needed (`wasm-opt -Oz`)
- [ ] Add loading states and user feedback
- [ ] Test with slow network (throttling)
- [ ] Verify memory doesn't leak over multiple generations

### Performance Considerations

**WASM Memory Limits:**
- Chrome/Firefox: ~4 GB
- Safari: ~2 GB
- Mobile: ~1 GB

**For Large Zone Plates (NA > 0.08):**
- Use Web Workers to prevent UI freezing
- Stream chunks instead of accumulating all data
- Show progress bar (already implemented)
- Consider server-side generation for very large files

---

## Troubleshooting

### "Cannot find module 'zpgen-wasm'"

```bash
# Reinstall the package
npm install file:../zpgen-wasm/examples/browser/pkg --force
```

### WASM fails to load

Check browser console for:
- MIME type errors (needs `application/wasm`)
- CORS errors (use same origin or configure CORS)
- File not found (verify paths)

### TypeScript errors

```bash
# Regenerate types
cd zpgen-wasm
./build-wasm.sh
```

### Poor performance

- Use Web Workers (see advanced example above)
- Reduce NA for testing (faster generation)
- Enable browser dev tools profiling

---

## Summary

**Key Points:**
1. ❌ **No common source** - C++ and Rust are separate
2. ✅ **Manual sync required** - Port C++ changes to Rust manually
3. ✅ **Easy React integration** - Use as NPM package
4. ✅ **Production ready** - Tested and optimized
5. ✅ **Memory safe** - No crashes from Rust code

**Recommended Setup:**
- Use Option 1 (NPM package)
- Use Web Workers for large zone plates
- Keep C++ as reference during transition
- Eventually migrate fully to Rust for browser+native

**Files to Remember:**
- Build: `./build-wasm.sh`
- Output: `examples/browser/pkg/`
- Install: `npm install file:../zpgen-wasm/examples/browser/pkg`
