import * as sgl from './sgl.js';

let canvas = document.getElementById("renderViewport");

let camera1 = new sgl.Camera(
    canvas,
    0.1, // near
    20, // far
    90, // fov
    new sgl.Vector3(0, 0, 0), // position
    new sgl.Vector3(0, 0, 0), // orientation
);

camera1.speed = 0.015;

let scene = new sgl.Scene(camera1);
let currentCamera = scene.currentCamera;

let renderer = new sgl.Renderer(canvas, scene);

let keys = {};
let gamepadIndex = null;

document.addEventListener('DOMContentLoaded', () => {
    canvas.addEventListener('click', () => {
        canvas.requestPointerLock();
        addEventListener('mousemove', handleMouseMove);
        addEventListener('keydown', handleKeyDown);
        addEventListener('keyup', handleKeyUp);
    });

    function handleMouseMove(event) {
        let mouseXMovement = event.movementX;
        let mouseYMovement = event.movementY;
    
        currentCamera.orientation.y -= mouseXMovement * currentCamera.speed * 10;
        currentCamera.orientation.x -= mouseYMovement * currentCamera.speed * 10;

        cameraClamp();
    }

    function handleKeyDown(event) {
        let key = event.key;

        if (key.length === 1) {
            key = key.toLowerCase();
        }
    
        keys[key] = true;
    }

    function handleKeyUp(event) {
        let key = event.key;

        if (key.length === 1) {
            key = key.toLowerCase();
        }
    
        keys[key] = false;
    }

    document.addEventListener('pointerlockchange', () => {
        if (document.pointerLockElement === null) {
            removeEventListener('mousemove', handleMouseMove);
            removeEventListener('keydown', handleKeyDown);
            removeEventListener('keyup', handleKeyUp);
        }
    });

    window.addEventListener("gamepadconnected", (event) => {
        console.log("Gamepad connected at index", event.gamepad.index);
        gamepadIndex = event.gamepad.index;
    });

    window.addEventListener("gamepaddisconnected", () => {
        console.log("Gamepad disconnected");
        gamepadIndex = null;
    });
});

function cameraClamp() {
    let distanceFrom90 = Math.abs(90 - currentCamera.orientation.x);
    let distanceFrom270 = Math.abs(270 - currentCamera.orientation.x);

    if (currentCamera.orientation.x > 90 && currentCamera.orientation.x < 270) {
        if (distanceFrom90 <= distanceFrom270) {
            currentCamera.orientation.x = 90;
        } else {
            currentCamera.orientation.x = 270;
        }
    }
}

function moveForward() {
    currentCamera.position.x -= currentCamera.speed * Math.sin(sgl.SGLMath.degToRad(currentCamera.orientation.y));
    currentCamera.position.z -= currentCamera.speed * Math.cos(sgl.SGLMath.degToRad(currentCamera.orientation.y));
}

function moveBackward() {
    currentCamera.position.x += currentCamera.speed * Math.sin(sgl.SGLMath.degToRad(currentCamera.orientation.y));
    currentCamera.position.z += currentCamera.speed * Math.cos(sgl.SGLMath.degToRad(currentCamera.orientation.y));
}

function moveLeft() {
    currentCamera.position.x -= currentCamera.speed * Math.sin(sgl.SGLMath.degToRad(currentCamera.orientation.y + 90));
    currentCamera.position.z += currentCamera.speed * Math.cos(sgl.SGLMath.degToRad(currentCamera.orientation.y - 90));
}

function moveRight() {
    currentCamera.position.x += currentCamera.speed * Math.sin(sgl.SGLMath.degToRad(currentCamera.orientation.y + 90));
    currentCamera.position.z -= currentCamera.speed * Math.cos(sgl.SGLMath.degToRad(currentCamera.orientation.y - 90));
}

function moveUp() {
    currentCamera.position.y += currentCamera.speed;
}

function moveDown() {
    currentCamera.position.y -= currentCamera.speed;
}

function handleInputs() {
    if (keys["w"]) {
        moveForward();
    }
    if (keys["s"]) {
        moveBackward();
    }
    if (keys["a"]) {
        moveLeft();
    }
    if (keys["d"]) {
        moveRight();
    }
    if (keys[" "]) {
        moveUp();
    }
    if (keys["Shift"]) {
        moveDown();
    }

    handleGamepadInputs();

    cameraOrientationValueReassignment();
}

function handleGamepadInputs() {
    if (gamepadIndex === null) return;

    let gamepad = navigator.getGamepads()[gamepadIndex];
    if (!gamepad) return;

    const DEAD_ZONE = 0.2; // Dead zone threshold

    // Move the camera with the left stick, applying the dead zone
    let leftStickX = gamepad.axes[0];
    let leftStickY = gamepad.axes[1];
    
    // Apply dead zone for left stick
    if (Math.abs(leftStickX) > DEAD_ZONE) {
        if (leftStickX < -DEAD_ZONE) moveLeft();
        if (leftStickX > DEAD_ZONE) moveRight();
    }
    if (Math.abs(leftStickY) > DEAD_ZONE) {
        if (leftStickY < -DEAD_ZONE) moveForward();
        if (leftStickY > DEAD_ZONE) moveBackward();
    }

    // Move up and down with trigger buttons
    if (gamepad.buttons[1].pressed) moveDown(); // B
    if (gamepad.buttons[0].pressed) moveUp();   // A

    // Rotate camera with the right stick, applying the dead zone
    let rightStickX = gamepad.axes[2];
    let rightStickY = gamepad.axes[3];

    // Apply dead zone for right stick
    if (Math.abs(rightStickX) > DEAD_ZONE) {
        currentCamera.orientation.y -= rightStickX * currentCamera.speed * 50;
    }
    if (Math.abs(rightStickY) > DEAD_ZONE) {
        currentCamera.orientation.x -= rightStickY * currentCamera.speed * 50;
    }

    cameraClamp();
}


function cameraOrientationValueReassignment() {
    currentCamera.orientation.x %= 360;
    currentCamera.orientation.y %= 360;
    currentCamera.orientation.z %= 360;

    if (currentCamera.orientation.x < 0) {
        currentCamera.orientation.x = 360 + currentCamera.orientation.x;
    }
    if (currentCamera.orientation.y < 0) {
        currentCamera.orientation.y = 360 + currentCamera.orientation.y;
    }
    if (currentCamera.orientation.z < 0) {
        currentCamera.orientation.z = 360 + currentCamera.orientation.z;
    }
}

function resizeCanvas() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
}

let lastFrameTime = performance.now();

let renderLoop = setInterval(function() {
    handleInputs();

    // Render
    scene.act();
    resizeCanvas();

    const currentFrameTime = performance.now();
    const deltaTime = (currentFrameTime - lastFrameTime) / 1000;
    //console.log("Render time: " + deltaTime + "ms");

    renderer.render();
}, 0);

// Assets
for (let i = 0; i < 1; i++) {
    scene.add(new sgl.TriangleMesh(new sgl.Vector3(0, 0, i + 2), new sgl.Vector3(1, 1, 1)));
}