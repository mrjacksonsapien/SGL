/**
 * This is a basic, barely-optimized 3D engine where code clarity has been prioritized over performance.
 * The goal of this library is to fully illustrate in a more human-readable manner all the steps of 3D rendering
 * pipeline by using only sequential instructions (CPU-based instructions) which means sacrificing the parallelism
 * of the GPU. Warning: This library is made for learning purpose of the 3D graphic pipeline and not for heavy intensive
 * scenes rendering. The CanvasRenderingContext2D API is used for the final render.
 *
 * Important note: You can use different data structures but most of the time, functions assume that you pass as an
 * argument a Float32Array and will directly modify the elements inside it.
 *
 * @author mrjacksonsapien
 */

/**
 * Math functions.
 */
export class SGLMath {
    /**
     * Converts degrees to radians.
     * @param degrees Unit in degrees.
     * @returns {number} Unit in radians.
     */
    static degToRad(degrees) {
        return degrees * (Math.PI / 180);
    }
    static cot(degrees) {
        return 1 / Math.tan(SGLMath.degToRad(degrees));
    }
}

/**
 * Takes care of rendering. Uses the Matrix class for vertices transformation.
 */
export class Renderer {
    /**
     * The size of one triangle inside the triangle data array.
     * @type {number}
     */
    static SIZEOF_TRIANGLE_DATA = 3;
    /**
     * The size of one vertex inside the vertex data array.
     * @type {number}
     */
    static SIZEOF_VERTEX_DATA = 4;

    /**
     * Creates an instance to render the scene on the canvas using CanvasRenderingContext2D API.
     * @param canvas HTML canvas.
     * @param scene SGL Scene object.
     */
    constructor(canvas, scene) {
        this.canvas = canvas;
        this.scene = scene;
    }

    /**
     * Converts the Scene object attached to the Renderer into 2 Float32Arrays stored in one JavaScript array.
     * @returns {Float32Array[]} Contains triangles data and vertices data. The triangles array contains the indices
     * to their corresponding vertices in the vertices data array.
     */
    convertSceneToFlatArrays() {
        let vertices = [];
        const triangles = [];
        const meshes = this.scene.meshes;

        // Get all the vertices from all the triangles in the scene
        for (let i = 0; i < meshes.length; i++) {
            const mesh = meshes[i];

            for (let j = 0; j < mesh.triangles.length; j++) {
                const triangle = mesh.triangles[j];
                triangles.push(triangle);

                vertices.push(triangle.vertex1);
                vertices.push(triangle.vertex2);
                vertices.push(triangle.vertex3);
            }
        }

        // Remove duplicates (many triangles sharing the same vertex)
        vertices = [...new Set(vertices)];

        const trianglesData = new Float32Array(triangles.length * Renderer.SIZEOF_TRIANGLE_DATA);

        // Map vertices to triangles
        for (let i = 0; i < triangles.length; i++) {
            const triangle = triangles[i];
            const triangleIndex = i * Renderer.SIZEOF_TRIANGLE_DATA;

            trianglesData[triangleIndex] = vertices.indexOf(triangle.vertex1) * Renderer.SIZEOF_VERTEX_DATA;
            trianglesData[triangleIndex + 1] = vertices.indexOf(triangle.vertex2) * Renderer.SIZEOF_VERTEX_DATA;
            trianglesData[triangleIndex + 2] = vertices.indexOf(triangle.vertex3) * Renderer.SIZEOF_VERTEX_DATA;
        }

        let verticesData = new Float32Array(vertices.length * Renderer.SIZEOF_VERTEX_DATA);

        // Add vertex values into flat array
        for (let i = 0; i < vertices.length; i++) {
            const vertex = vertices[i];
            const verticesDataIndex = i * Renderer.SIZEOF_VERTEX_DATA;

            verticesData[verticesDataIndex] = vertex.position.x;
            verticesData[verticesDataIndex + 1] = vertex.position.y;
            verticesData[verticesDataIndex + 2] = vertex.position.z;
            verticesData[verticesDataIndex + 3] = 1;
        }

        return [trianglesData, verticesData];
    }

    /**
     * Apply a matrix to a flat array of many vertices (Converts each vertex individually).
     * @param matrix The matrix obtained from a method in the Matrix class.
     * @param verticesData The vertices data float array.
     */
    applyMatrixToVertices(matrix, verticesData) {
        for (let i = 0; i < verticesData.length / Renderer.SIZEOF_VERTEX_DATA; i++) {
            Matrix.multiplyMatrixWithVertex(matrix,  i * Renderer.SIZEOF_VERTEX_DATA, verticesData);
        }
    }

    /**
     * Clips all the triangles based on 3 situations and returns a new JavaScript array containing new data for the
     * triangles and vertices (same structure as return of convertSceneToFlatArrays()).
     * A: Completely outside. Discards the triangle.
     * B: Completely inside. Keeps the triangle as-is.
     * C: Some vertices are outside a plane. clip to the plane.
     * @param trianglesData The triangles data.
     * @param verticesData The vertices data.
     * @returns {Float32Array[]} The new JavaScript array containing the triangles and vertices data.
     */
    clip(trianglesData, verticesData) {
        const newTriangles = [];
        const newVertices = [];

        const indexMap = {};

        function addVertex(vertexIndex) {
            if (indexMap[vertexIndex] === undefined) {
                indexMap[vertexIndex] = newVertices.length;

                newVertices.push(
                    verticesData[vertexIndex],
                    verticesData[vertexIndex + 1],
                    verticesData[vertexIndex + 2],
                    verticesData[vertexIndex + 3]
                );
            }
            return indexMap[vertexIndex];
        }

        function clipTriangle(triangleIndex) {
            function planes(vertexIndex) {
                const x = verticesData[vertexIndex];
                const y = verticesData[vertexIndex + 1];
                const z = verticesData[vertexIndex + 2];
                const w = verticesData[vertexIndex + 3];

                return [
                    x >= w, // Left plane
                    x <= -w, // Right plane
                    y >= w, // Top plane
                    y <= -w, // Bottom plane
                    z >= w, // Near plane
                    z <= -w // Far plane
                ]
            }

            function isInside(vertexPlanes) {
                return vertexPlanes.every(flag => flag);
            }


            function clipAgainstPlanes(v1, v2, v3, v1Planes, v2Planes, v3Planes) {
                for (let i = 0; i < 6; i++) {
                    const v1IsInside = v1Planes[i];
                    const v2IsInside = v2Planes[i];
                    const v3IsInside = v3Planes[i];

                    const insideVertices = [];
                    const outsideVertices = [];

                    if (v1IsInside) insideVertices.push(v1);
                    else outsideVertices.push(v1);

                    if (v2IsInside) insideVertices.push(v2);
                    else outsideVertices.push(v2);

                    if (v3IsInside) insideVertices.push(v3);
                    else outsideVertices.push(v3);

                    if (outsideVertices.length === 2) {
                        // TODO: Two vertices outside
                    } else {
                        // TODO: One vertex outside
                    }
                }
            }

            const v1 = trianglesData[triangleIndex];
            const v2 = trianglesData[triangleIndex + 1];
            const v3 = trianglesData[triangleIndex + 2];

            const v1Planes = planes(v1);
            const v2Planes = planes(v2);
            const v3Planes = planes(v3);

            // Keep triangle as-is if all vertices are inside
            if (isInside(v1Planes) && isInside(v2Planes) && isInside(v3Planes)) {
                newTriangles.push(
                    addVertex(v1),
                    addVertex(v2),
                    addVertex(v3)
                );
                return;
            }

            // Discard triangle if all vertices are outside
            for (let i = 0; i < 6; i++) {
                if (!v1Planes[i] && !v2Planes[i] && !v3Planes[i]) {
                    return;
                }
            }

            // Clip triangle if one or two vertices are outside
            clipAgainstPlanes(v1, v2, v3, v1Planes, v2Planes, v3Planes);
        }

        for (let i = 0; i < trianglesData.length / Renderer.SIZEOF_TRIANGLE_DATA; i++) {
            clipTriangle(i * Renderer.SIZEOF_TRIANGLE_DATA);
        }

        return [new Float32Array(newTriangles), new Float32Array(newVertices)];
    }

    /**
     * Applies perspective division to clip-space vertices.
     * @param verticesData Clip-space vertices.
     */
    applyPerspectiveDivisionToClipVertices(verticesData) {
        for (let i = 0; i < verticesData.length / Renderer.SIZEOF_VERTEX_DATA; i++) {
            const vertexIndex = i * Renderer.SIZEOF_VERTEX_DATA;

            verticesData[vertexIndex] /= verticesData[vertexIndex + 3];
            verticesData[vertexIndex + 1] /= verticesData[vertexIndex + 3];
            verticesData[vertexIndex + 2] /= verticesData[vertexIndex + 3];
        }
    }

    /**
     * Maps the NDC vertices to the canvas 2D coordinates.
     * @param verticesData NDC vertices.
     */
    mapNdcVerticesToScreenCoordinates(verticesData) {
        for (let i = 0; i < verticesData.length / Renderer.SIZEOF_VERTEX_DATA; i++) {
            const vertexIndex = i * Renderer.SIZEOF_VERTEX_DATA;

            verticesData[vertexIndex] = (verticesData[vertexIndex] + 1) / 2 * this.canvas.clientWidth;
            verticesData[vertexIndex + 1] = (1 - verticesData[vertexIndex + 1]) / 2 * this.canvas.clientHeight;
        }
    }

    /**
     * Display the triangle's wireframe.
     * @param trianglesData Triangles data after clipping.
     * @param verticesData Vertices data after 2D mapping.
     */
    displayTrianglesWireframe(trianglesData, verticesData) {
        const ctx = this.canvas.getContext('2d');
        for (let i = 0; i < trianglesData.length / Renderer.SIZEOF_TRIANGLE_DATA; i++) {
            const triangleIndex = i * Renderer.SIZEOF_TRIANGLE_DATA;

            const vertex1Index = trianglesData[triangleIndex];
            const vertex2Index = trianglesData[triangleIndex + 1];
            const vertex3Index = trianglesData[triangleIndex + 2];

            ctx.beginPath();
            ctx.moveTo(verticesData[vertex1Index], verticesData[vertex1Index + 1]);
            ctx.lineTo(verticesData[vertex2Index], verticesData[vertex2Index + 1]);
            ctx.lineTo(verticesData[vertex3Index], verticesData[vertex3Index + 1]);
            ctx.closePath();
            ctx.strokeStyle = 'black';
            ctx.stroke();
        }
    }

    /**
     * Goes through all the steps of the 3D rendering pipeline. The scene must have a current camera when this method is called.
     * Since JavaScript is really easy to break, make sure that all your attributes are initialized and defined. If you use the classes
     * correctly, it shouldn't be a problem but be careful.
     */
    render() {
        const camera = this.scene.currentCamera;

        let sceneData = this.convertSceneToFlatArrays();
        let trianglesData = sceneData[0];
        let verticesData = sceneData[1];

        /* Now in world space (coordinates relative to the world origin). The scene has been converted into
        2 flat (one dimension) arrays: One for all the vertices in the scene, and one for the triangles. The triangles
        hold the indices (in verticesData) of the vertices that shapes them. You can see this as some kind of "reference".
        Making flat arrays will help make calculations faster than accessing objects data. */

        this.applyMatrixToVertices(Matrix.createViewMatrix(camera), verticesData);

        // Now in view space: Peform calculations for lighting, normal vectors, etc.

        this.applyMatrixToVertices(Matrix.createProjectionMatrix(camera), verticesData);

        /* Now in clip space: Perform clipping and triangles subdivision. The variables holding the data arrays are updated
        because some triangles and vertices might have been discarded with some new ones added, which requires changing the size
        of the arrays. */

        sceneData = this.clip(trianglesData, verticesData);
        trianglesData = sceneData[0];
        verticesData = sceneData[1];

        this.applyPerspectiveDivisionToClipVertices(verticesData);

        /* Now in NDC space. The perspective division is the part where everything has been turned based on the perspective
        effect (things further away are smaller/vertices are closer to the center of the viewport) */

        this.mapNdcVerticesToScreenCoordinates(verticesData);

        /* Now in 2D screen space. The vertices are finally mapped to the 2D screen coordinates. In this space, we can finally
        see the final position of the vertices. */

        this.displayTrianglesWireframe(trianglesData, verticesData);

        console.log(verticesData.length / Renderer.SIZEOF_VERTEX_DATA, trianglesData.length / Renderer.SIZEOF_TRIANGLE_DATA);
    }
}

/**
 * Holds all the hardcoded matrices for 3D transformation. Some conventional aspects are respected while others like
 * matrix multiplications might change from one engine to another.
 */
export class Matrix {
    /**
     * Multiplies a 4x4 matrix (flat array) with a vertex.
     * @param m The matrix (flat array).
     * @param vertexIndex The index of the first component of the vertex (x) inside the verticesData.
     * @param verticesData The vertices data (Float32Array).
     */
    static multiplyMatrixWithVertex(m, vertexIndex, verticesData) {
        const x = verticesData[vertexIndex];
        const y = verticesData[vertexIndex + 1];
        const z = verticesData[vertexIndex + 2];
        const w = verticesData[vertexIndex + 3];

        verticesData[vertexIndex] = x * m[0] + y * m[4] + z * m[8] + w * m[12];
        verticesData[vertexIndex + 1] = x * m[1] + y * m[5] + z * m[9] + w * m[13];
        verticesData[vertexIndex + 2] = x * m[2] + y * m[6] + z * m[10] + w * m[14];
        verticesData[vertexIndex + 3] = x * m[3] + y * m[7] + z * m[11] + w * m[15];
    }

    /**
     * Multiplies two 4x4 matrices together.
     * @param a First matrix (flat array).
     * @param b Second matrix (flat array).
     * @returns {Float32Array} Result 4x4 matrix.
     */
    static multiply4x4Matrices(a, b) {
        const result = new Float32Array(16);

        for (let i = 0; i < 4; i++) {
            const offsetA = i * 4;
            for (let j = 0; j < 4; j++) {
                const offsetB = j;
                result[offsetA + j] = a[offsetA] * b[offsetB] +
                    a[offsetA + 1] * b[4 + offsetB] +
                    a[offsetA + 2] * b[8 + offsetB] +
                    a[offsetA + 3] * b[12 + offsetB];
            }
        }

        return result;
    }

    /**
     * Creates a translation matrix with x, y, z position components.
     * @param x
     * @param y
     * @param z
     * @returns {Float32Array}
     */
    static createTranslationMatrix(x, y, z) {
        const translationMatrix = new Float32Array(16);
        translationMatrix.set([
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            x, y, z, 1
        ]);
        return translationMatrix;
    }

    /**
     * Creates a pitch matrix (x rotation).
     * @param pitchDegrees Degrees for x.
     * @returns {Float32Array}
     */
    static createPitchMatrix(pitchDegrees) {
        const pitch = SGLMath.degToRad(pitchDegrees);
        const pitchMatrix = new Float32Array(16);
        pitchMatrix.set([
            1, 0, 0, 0,
            0, Math.cos(pitch), -Math.sin(pitch), 0,
            0, Math.sin(pitch), Math.cos(pitch), 0,
            0, 0, 0, 1
        ]);
        return pitchMatrix;
    }

    /**
     * Create a yaw matrix (y).
     * @param yawDegrees Degrees for y.
     * @returns {Float32Array}
     */
    static createYawMatrix(yawDegrees) {
        const yaw = SGLMath.degToRad(yawDegrees);
        const yawMatrix = new Float32Array(16);
        yawMatrix.set([
            Math.cos(yaw), 0, Math.sin(yaw), 0,
            0, 1, 0, 0,
            -Math.sin(yaw), 0, Math.cos(yaw), 0,
            0, 0, 0, 1
        ]);
        return yawMatrix;
    }

    /**
     * Creates a roll matrix (z).
     * @param rollDegrees Degrees for z.
     * @returns {Float32Array}
     */
    static createRollMatrix(rollDegrees) {
        const roll = SGLMath.degToRad(rollDegrees);
        const rollMatrix = new Float32Array(16);
        rollMatrix.set([
            Math.cos(roll), -Math.sin(roll), 0, 0,
            Math.sin(roll), Math.cos(roll), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        ]);
        return rollMatrix;
    }

    /**
     * Creates an euleur matrix respecting the video game camera rotation standard.
     * @param x In degrees.
     * @param y In degrees.
     * @param z In degrees.
     * @returns {Float32Array}
     */
    static createEulerMatrix(x, y, z) {
        return Matrix.multiply4x4Matrices(
            Matrix.multiply4x4Matrices(
                Matrix.createYawMatrix(y),
                Matrix.createPitchMatrix(x)
            ),
            Matrix.createRollMatrix(z)
        );
    }

    /**
     * Creates a view matrix for the camera.
     * @param camera Camera object.
     * @returns {Float32Array}
     */
    static createViewMatrix(camera) {
        return Matrix.multiply4x4Matrices(
            Matrix.createTranslationMatrix(camera.position.x, camera.position.y, camera.position.z),
            Matrix.createEulerMatrix(camera.orientation.x, camera.orientation.y, camera.orientation.z)
        )
    }

    /**
     * Creates a projection matrix (clip-space transformation).
     * @param camera Camera object.
     * @returns {Float32Array}
     */
    static createProjectionMatrix(camera) {
        const aspectRatio = camera.canvas.clientWidth / camera.canvas.clientHeight;
        const projectionMatrix = new Float32Array(16);
        projectionMatrix.set([
            SGLMath.cot(camera.fov / 2) / aspectRatio, 0, 0, 0,
            0, SGLMath.cot(camera.fov / 2), 0, 0,
            0, 0, -(camera.far / (camera.far - camera.near)), -1,
            0, 0, (camera.far * camera.near) / (camera.far - camera.near), 0
        ]);
        return projectionMatrix;
    }
}

/**
 * Contains all the meshes that will be rendered and the current camera that will be used to render the scene.
 */
export class Scene {
    /**
     * The camera that will be set as the current camera for rendering.
     * @param camera Camera object.
     */
    constructor(camera) {
        this.currentCamera = camera;
        this.meshes = [];
    }

    /**
     * Add a mesh to the scene.
     * @param mesh
     */
    add(mesh) {
        this.meshes.push(mesh);
    }

    /**
     * Removes a mesh from the scene.
     * @param mesh
     */
    remove(mesh) {
        let index = this.meshes.indexOf(mesh);
        if (index !== -1) {
            this.meshes.splice(index, 1);
        }
    }

    /**
     * Logic to update the scene.
     */
    act() {
        this.meshes.forEach(mesh => {
            mesh.act();
        });
    }
}

/**
 * Camera object used for rendering a scene.
 */
export class Camera {
    /**
     * Creates the camera.
     * @param canvas HTML canvas used for rendering the scene.
     * @param near The near clipping plane (cannot be equal or under 0).
     * @param far Has to be larger than the near value.
     * @param fov Field of view.
     * @param position Position of the camera.
     * @param orientation Orientation of the camera.
     */
    constructor(canvas, near, far, fov, position, orientation) {
        this.canvas = canvas;
        this.near = near;
        this.far = far;
        this.fov = fov;
        this.position = position;
        this.orientation = orientation;
    }
}

/**
 * Mesh object containing all the data for rendering.
 */
export class Mesh {
    constructor(triangles) {
        this.triangles = triangles
        this.act = () => {};
    }
}

/**
 * Contains the references to 3 Vertex and the data of the triangle.
 */
export class Triangle {
    constructor(vertex1, vertex2, vertex3) {
        this.vertex1 = vertex1;
        this.vertex2 = vertex2;
        this.vertex3 = vertex3;
    }
}

/**
 * Contains the data of a vertex.
 */
export class Vertex {
    constructor(position) {
        this.position = position;
    }
}

/**
 * Contains an x, y and z component.
 */
export class Vector3 {
    constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
}

// Templates //

export class Cube extends Mesh {
    constructor(position, size) {
        super();
        this.position = position;
        this.size = size;

        let vertices = [
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y + this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y + this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y - this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y - this.size.y / 2, this.position.z - this.size.z / 2)),
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y + this.size.y / 2, this.position.z + this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y + this.size.y / 2, this.position.z + this.size.z / 2)),
            new Vertex(new Vector3(this.position.x + this.size.x / 2, this.position.y - this.size.y / 2, this.position.z + this.size.z / 2)),
            new Vertex(new Vector3(this.position.x - this.size.x / 2, this.position.y - this.size.y / 2, this.position.z + this.size.z / 2))
        ];

        this.triangles = [
            new Triangle(vertices[0], vertices[1], vertices[2]),
            new Triangle(vertices[0], vertices[2], vertices[3]),
            new Triangle(vertices[4], vertices[5], vertices[6]),
            new Triangle(vertices[4], vertices[6], vertices[7]),
            new Triangle(vertices[0], vertices[1], vertices[5]),
            new Triangle(vertices[0], vertices[5], vertices[4]),
            new Triangle(vertices[3], vertices[2], vertices[6]),
            new Triangle(vertices[3], vertices[6], vertices[7]),
            new Triangle(vertices[0], vertices[4], vertices[3]),
            new Triangle(vertices[4], vertices[3], vertices[7]),
            new Triangle(vertices[1], vertices[5], vertices[6]),
            new Triangle(vertices[1], vertices[2], vertices[6])
        ];
    }
}