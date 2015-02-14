/*
    Copyright (C) 2014-2015 Nate Craun

    This file is part of Jubal.

    Jubal is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jubal is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jubal.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "test.h"
#include "lfqueue.h"

int main(int argc, char **argv) {

    {
        Lfqueue *q = make_lfqueue(10, sizeof(int));

        TEST("constructor", lfqueue_is_empty(q));
        TEST("constructor", q->capacity == 10);
        free_lfqueue(q);
    }

    {
        Lfqueue *q = make_lfqueue(6, sizeof(int));
        for (int i = 0; i < 5; i++) {
            TEST("normal add remove", lfqueue_size(q) == i);
            lfqueue_add(q, &i);
            TEST("normal add remove", lfqueue_size(q) == i+1);
        }

        for (int i = 0; i < 5; i++) {
            TEST("normal add remove", lfqueue_size(q) == 5-i);
            int out;
            lfqueue_poll(q, &out);
            TEST("normal add remove", out == i);
            TEST("normal add remove", lfqueue_size(q) == 5-i-1);
        }
        free_lfqueue(q);
    }

    {
        Lfqueue *q = make_lfqueue(6, sizeof(int));
        for (int i = 0; i < 5; i++) {
            lfqueue_add(q, &i);
            TEST("interleaved add remove", lfqueue_size(q) == 1);
            lfqueue_advance(q);
            TEST("interleaved add remove", lfqueue_size(q) == 0);
        }
        TEST("interleaved add remove", lfqueue_size(q) == 0);
        TEST("interleaved add remove", q->capacity == 6);
        free_lfqueue(q);
    }

    {
        Lfqueue *q = make_lfqueue(6,sizeof(int));
        int i0 = 0;
        lfqueue_add(q, &i0);
        int i1 = 1;
        lfqueue_add(q, &i1);

        TEST("peek", lfqueue_size(q) == 2);

        for (int i = 0; i < 3; i++) {
            int out;
            lfqueue_peek(q, &out);
            TEST("peek", out == 0);
            TEST("peek", lfqueue_size(q) == 2);
        }
        free_lfqueue(q);
    }

    {
        Lfqueue *q = make_lfqueue(6, sizeof(int));
/*
        {
            Option<int> o = q.pop_front();
            TEST("pop from empty", !o);
        }

        {
            Option<int> o = q.peek_front();
            TEST("peek from empty", !o);
        }
*/
        {
            for (int i = 0; i < 5; i++) {
                int i0 = 0;
                lfqueue_add(q, &i0);
            }
/*
            TEST("push to full", lfqueue_size(q) == 5);
            TEST("push to full", lfqueue_size(q) == q->capacity);

            bool push_succeed = lfqueue_add(q, 1);
            TEST("push to full", !push_succeed);
            TEST("push to full", lfqueue_size(q) == 5);
            TEST("push to full", lfqueue_size(q) == q->capacity);
*/
        }
        free_lfqueue(q);
    }

   return 0; 
}
